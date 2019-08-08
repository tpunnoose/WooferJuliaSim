using OSQP
# using Parametron # TODO: use this

@with_kw struct MPCControllerParams
	# gravity
	g::Float64 = 9.81

	min_vert_force::Float64 = 1.0
	max_vert_force::Float64 = 133.0

	# inverse of inertia tensor in body frame
	J_inv::Diagonal{Float64}

	# discretization length
	dt::Float64

	# planning horizon length
	N::Int64

	A_c::SparseMatrixCSC{Float64, Int64}

	# system matrices
	A_qp::Matrix{Float64} = zeros(N*10, 10)
	B_qp::Matrix{Float64} = zeros(N*10, N*12)

	cont_sys::Matrix{Float64} = zeros(22, 22)
	disc_sys::Matrix{Float64} = zeros(22, 22)
	B_c::Matrix{Float64} = zeros(10, 12)
	r_hat::Matrix{Float64} = zeros(3,3)

	A_d::Matrix{Float64} = zeros(10, 10)
	B_d_vec::Array{Array{Float64,2}, 1} = [zeros(10,12) for i in 1:N]

	C::SparseMatrixCSC{Float64, Int64}
	lb::Vector{Float64}
	ub::Vector{Float64}

	K::Diagonal{Float64}
	R::Diagonal{Float64}
	alpha::Float64 = 1e-6

	# OSQP model
	prob::OSQP.Model
	first_run::Vector{Bool} = [true]
	P::SparseMatrixCSC{Float64, Int64} = zeros(12*N, 12*N)
	q::Vector{Float64} = zeros(12*N)
end

function initMPCControllerConfig(dt::Number, N::Integer, woofer_config::WooferConfig)
	J_inv = inv(woofer_config.INERTIA)

	A_c = zeros(10,10)
	A_c[1:3, 6:8] = Matrix{Float64}(I, 3, 3)
	A_c[6, 10] = -1
	A_c = sparse(A_c)

	prob = OSQP.Model()

	# TODO: construct constant constraint matrices here
	C_i = zeros(20, 12)
	mu = 0.6
	min_vert_force = 1
	max_vert_force = 133

	for i in 1:4
		# fz >= 0
		C_i[(i-1)*5+1, (i-1)*3+3] = 1
		# ufz+fx >= 0
		C_i[(i-1)*5+2, (i-1)*3+1] = 1
		C_i[(i-1)*5+2, (i-1)*3+3] = mu
		# fx-ufz <= 0
		C_i[(i-1)*5+3, (i-1)*3+1] = 1
		C_i[(i-1)*5+3, (i-1)*3+3] = -mu
		# ufz+fy >= 0
		C_i[(i-1)*5+4, (i-1)*3+2] = 1
		C_i[(i-1)*5+4, (i-1)*3+3] = mu
		# fy-ufz <= 0
		C_i[(i-1)*5+5, (i-1)*3+2] = 1
		C_i[(i-1)*5+5, (i-1)*3+3] = -mu
	end

	lb_i = zeros(20)
	ub_i = zeros(20)

	for i in 1:4
		# fz >= 0
		lb_i[(i-1)*5+1] = min_vert_force
		ub_i[(i-1)*5+1] = max_vert_force
		# ufz+fx >= 0
		lb_i[(i-1)*5+2] = 0
		ub_i[(i-1)*5+2] = Inf
		# fx-ufz <= 0
		lb_i[(i-1)*5+3] = -Inf
		ub_i[(i-1)*5+3] = 0
		# ufz+fy >= 0
		lb_i[(i-1)*5+4] = 0
		ub_i[(i-1)*5+4] = Inf
		# fy-ufz >= 0
		lb_i[(i-1)*5+5] = -Inf
		ub_i[(i-1)*5+5] = 0
	end

	C = zeros(20*N, 12*N)
	for i in 1:N
		C[20*(i-1)+1:20*(i-1)+20, 12*(i-1)+1:12*(i-1)+12] .= C_i
	end
	C = sparse(C)

	lb = repeat(lb_i, N)
	ub = repeat(ub_i, N)

	k = [50, 10, 50, 1, 1, 1, 1, 1, 1, 0]
	K = Diagonal(repeat(k, N))
	R = Diagonal(repeat(ones(12), N))

	mpcConfig = MPCControllerParams(prob=prob, J_inv=J_inv, dt=dt, N=N, K=K, R=R, A_c=A_c, C=C, lb=lb, ub=ub, min_vert_force=min_vert_force, max_vert_force=max_vert_force)
	return mpcConfig
end

function generateReferenceTrajectory!(x_ref::Array{T, 2}, x_curr::Vector{T}, x_des::Vector{T}, mpc_config::MPCControllerParams) where {T<:Number}
	# how agressively the MPC will attempt to get to desired state
	stiffness = 1.0

	N_ = floor(stiffness*mpc_config.N)

	x_diff = 1/N_ * (x_des - x_curr)
	for i in 1:mpc_config.N
		if i<N_
			x_ref[:, i] .= x_curr + i*x_diff
		else
			x_ref[:, i] .= x_des
		end
	end
end

function skewSymmetricMatrix!(A::Matrix, a::Vector)
	A[1,2] = -a[3]
	A[1,3] = a[2]
	A[2,1] = a[3]
	A[2,3] = -a[1]
	A[3,1] = -a[2]
	A[3,2] = a[1]
end

function solveFootForces!(forces::Vector{T}, x0::Vector{T}, x_ref::Array{T, 2}, contacts::Array{Int,2}, foot_locs::Array{T,2}, mpc_config::MPCControllerParams, woofer_config::WooferConfig) where {T<:Number}
	# x_ref: 10xN matrix of state reference trajectory
	# contacts: 4xN matrix of foot contacts over the planning horizon
	# foot_locs: 12xN matrix of foot location in body frame over planning horizon

	for i in 1:mpc_config.N
		## construct continuous time B matrix
		for j in 1:4
			if contacts[j,i] == 1
				mpc_config.B_c[4:6, (3*(j-1)+1):(3*(j-1)+3)] .= 1/woofer_config.MASS*contacts[j,i]*Matrix{Float64}(I, 3, 3)

				skewSymmetricMatrix!(mpc_config.r_hat, foot_locs[(3*(j-1)+1):(3*(j-1)+3),i])
				mpc_config.B_c[7:9, (3*(j-1)+1):(3*(j-1)+3)] .= mpc_config.J_inv*mpc_config.r_hat

				mpc_config.lb[20*(i-1)+ 1 + (j-1)*5+1] = mpc_config.min_vert_force
				mpc_config.ub[20*(i-1)+ 1 + (j-1)*5+1] = mpc_config.max_vert_force
			else
				mpc_config.B_c[4:9, (3*(j-1)+1):(3*(j-1)+3)] .= zeros(6, 3)

				mpc_config.lb[20*(i-1)+ 1 + (j-1)*5+1] = 0
				mpc_config.ub[20*(i-1)+ 1 + (j-1)*5+1] = 0
			end
		end

		## get the discretized system matrices via ZOH approximation:
		mpc_config.cont_sys[1:10, 1:10] .= mpc_config.A_c
		mpc_config.cont_sys[1:10, 11:22] .= mpc_config.B_c

		mpc_config.disc_sys .= exp(mpc_config.cont_sys*mpc_config.dt) # TODO: see if there is a better way to do this
		if i == 1
			mpc_config.A_d .= mpc_config.disc_sys[1:10, 1:10]
		end
		mpc_config.B_d_vec[i] .= mpc_config.disc_sys[1:10, 11:22]

		## add row of A_qp
		mpc_config.A_qp[(10*(i-1)+1):(10*(i-1)+10), :] .= mpc_config.A_d^i

		## add row of B_qp
		for l in 0:(i-1)
			mpc_config.B_qp[(10*(i-1)+1):(10*(i-1)+10), (12*l+1):(12*l+12)] .= mpc_config.A_d^(i-1-l)*mpc_config.B_d_vec[l+1]
		end
	end

	mpc_config.P .= sparse(2*(mpc_config.B_qp' * mpc_config.K * mpc_config.B_qp + mpc_config.alpha*mpc_config.R))
	mpc_config.q .= sparse(2*(mpc_config.B_qp' * mpc_config.K * mpc_config.A_qp * x0 - mpc_config.B_qp' * mpc_config.K * reshape(x_ref, 10*mpc_config.N)))

	# TODO: make this work
	# if mpc_config.first_run[1]
	# 	OSQP.setup!(mpc_config.prob; P=mpc_config.P, q=mpc_config.q, A=mpc_config.C, l=mpc_config.lb, u=mpc_config.ub, verbose=false)
	# 	mpc_config.first_run[1] = false
	# else
	# 	OSQP.update!(mpc_config.prob; q=mpc_config.q)
	# 	OSQP.update!(mpc_config.prob; Px=mpc_config.P)
	# end

	OSQP.setup!(mpc_config.prob; P=mpc_config.P, q=mpc_config.q, A=mpc_config.C, l=mpc_config.lb, u=mpc_config.ub, verbose=false)
	results = OSQP.solve!(mpc_config.prob)

	forces .= (results.x)[1:12]

	x_exp = mpc_config.A_qp*x0 + mpc_config.B_qp*results.x

	# println(x_exp[10*(mpc_config.N-1)+1:10*(mpc_config.N-1)+10])
end
