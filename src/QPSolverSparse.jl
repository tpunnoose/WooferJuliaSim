"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

using OSQP
using ControlSystems
# using Parametron # TODO: use this

@with_kw struct MPCControllerParams
	# inverse of inertia tensor in body frame
	J_inv::Diagonal{Float64}

	# discretization length
	dt::Float64

	# planning horizon length
	N::Int64

	B_c::Matrix{Float64} = zeros(12, 12)
	r_hat::Matrix{Float64} = zeros(3,3)

	A_d::Matrix{Float64} = zeros(12, 12)
	B_d::Matrix{Float64} = zeros(12, 12)

	C_dense::Matrix{Float64}
	# C_sparse::SparseMatrixCSC{Float64, Int64} # TODO: Add this
	lb::Vector{Float64}
	ub::Vector{Float64}

	Q::Diagonal{Float64}
	R::Diagonal{Float64}
	H::Array{Float64, 2}
	V::Array{Float64, 2} = zeros(12,12)

	# sections of constraint matrix C
	B_lineq::Matrix{Float64}

	u_ref::Vector{Float64}
	z_ref::Vector{Float64} = zeros(24*N + 12)
	z_result::Vector{Float64} = zeros(24*N + 12)

	# OSQP model
	prob::OSQP.Model
	first_run::Vector{Bool} = [true]
	P::SparseMatrixCSC{Float64, Int64} = zeros(24*N+12, 24*N+12)
	q::Vector{Float64} = zeros(24*N+12)
end

function initMPCControllerConfig(dt::Number, N::Integer)
	# z = (x0, ..., xN, u0, ..., u0, ..., uN)
	J_inv = inv(WOOFER_CONFIG.INERTIA)

	A_c = zeros(12,12)
	A_c[1:3, 7:9] = Matrix{Float64}(I, 3, 3)
	A_c[4:6, 10:12] = Matrix{Float64}(I, 3, 3)

	prob = OSQP.Model()

	# TODO: construct constant constraint matrices here
	C_i = zeros(20, 12)
	mu = 1.5
	min_vert_force = 1
	max_vert_force = 133

	# friction constraint matrix for each u_i
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

	lb_friction = zeros(20)
	ub_friction = zeros(20)

	# friction lower and upper bounds for each u_i
	for i in 1:4
		# # fz >= 0
		# lb_friction[(i-1)*5+1] = min_vert_force
		# ub_friction[(i-1)*5+1] = max_vert_force
		# # ufz+fx >= 0
		# lb_friction[(i-1)*5+2] = 0
		# ub_friction[(i-1)*5+2] = Inf
		# # fx-ufz <= 0
		# lb_friction[(i-1)*5+3] = -Inf
		# ub_friction[(i-1)*5+3] = 0
		# # ufz+fy >= 0
		# lb_friction[(i-1)*5+4] = 0
		# ub_friction[(i-1)*5+4] = Inf
		# # fy-ufz >= 0
		# lb_friction[(i-1)*5+5] = -Inf
		# ub_friction[(i-1)*5+5] = 0
		# fz >= 0
		lb_friction[(i-1)*5+1] = -Inf
		ub_friction[(i-1)*5+1] = Inf
		# ufz+fx >= 0
		lb_friction[(i-1)*5+2] = -Inf
		ub_friction[(i-1)*5+2] = Inf
		# fx-ufz <= 0
		lb_friction[(i-1)*5+3] = -Inf
		ub_friction[(i-1)*5+3] = Inf
		# ufz+fy >= 0
		lb_friction[(i-1)*5+4] = -Inf
		ub_friction[(i-1)*5+4] = Inf
		# fy-ufz >= 0
		lb_friction[(i-1)*5+5] = -Inf
		ub_friction[(i-1)*5+5] = Inf
	end

	# gravity affine term -> maps to vdot
	d = [0, 0, 0, 0, 0, 0, 0, 0, -9.81, 0, 0, 0]*dt

	C_dense = zeros((20 + 12)*N + 12, 24*N+12)

	# constraint matrix for discrete LDS (enforces Ax[k-1] - x[k] + Bu[k-1] = 0)
	A_d = exp(A_c*dt)
	A_lineq = kron(Matrix{Float64}(I, N+1, N+1), -Matrix{Float64}(I, 12, 12)) + kron(diagm(-1 => ones(N)), A_d)

	# vector of B_d changes every horizon
	B_lineq = zeros(12*N, 12*N)

	# constant friction constraint matrix
	C_friction = kron(Matrix{Float64}(I, N, N), C_i)

	# OSQP Constraint bounds
	lb = vcat(zeros(12), kron(ones(N), -d), kron(ones(N), lb_friction))
	ub = vcat(zeros(12), kron(ones(N), -d), kron(ones(N), ub_friction))


	C_dense = zeros(32*N+12, 24*N+12)
	C_dense[1:12*(N+1), 1:(12*N+12)] .= A_lineq
	C_dense[(12*(N+1)+1):(32*N+12), (12*N+13):(24*N+12)] .= C_friction

	q = [1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e3, 1e3, 1e3, 1, 1, 1e2]
	r = [1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4]

	Q = Diagonal(q)
	R = Diagonal(r)

	H = Diagonal(vcat(repeat(q, N+1), repeat(r, N)))

	u_ref = [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4

	mpcConfig = MPCControllerParams(J_inv=J_inv, dt=dt, N=N, H=H, A_d=A_d, C_dense=C_dense, B_lineq=B_lineq, lb=lb, ub=ub, prob=prob, u_ref=u_ref, Q=Q, R=R)
	return mpcConfig
end

function generateReferenceTrajectory!(x_ref::Array{T, 2}, x_curr::Vector{T}, x_des::Vector{T}, mpc_config::MPCControllerParams) where {T<:Number}
	# TODO: integrate the x,y,ψ position from the reference

	x_diff = 1/mpc_config.N * (x_des - x_curr)
	for i in 1:mpc_config.N
		if i<mpc_config.N
			x_ref[:, i] .= x_curr + i*x_diff
		else
			x_ref[:, i] .= x_des
		end
	end
end

function solveFootForces!(forces::Vector{T}, x0::Vector{T}, x_ref::Array{T, 2}, contacts::Array{Int,2}, foot_locs::Array{T,2}, mpc_config::MPCControllerParams, woofer_config::WooferConfig) where {T<:Number}
	# x_ref: 12xN matrix of state reference trajectory
	# contacts: 4xN matrix of foot contacts over the planning horizon
	# foot_locs: 12xN matrix of foot location in body frame over planning horizon

	# average yaw over the reference trajectory -> TODO: you don't have to do this in sparse formulation
	ψ = sum(x_ref[6,:])/mpc_config.N
	# TODO: recalculate A_c(ψ) here

	for i in 1:mpc_config.N
		## construct continuous time B matrix
		for j in 1:4
			if contacts[j,i] == 1
				mpc_config.B_c[7:9, (3*(j-1)+1):(3*(j-1)+3)] .= 1/WOOFER_CONFIG.MASS*Matrix{Float64}(I, 3, 3)

				skewSymmetricMatrix!(mpc_config.r_hat, foot_locs[(3*(j-1)+1):(3*(j-1)+3)])
				mpc_config.B_c[10:12, (3*(j-1)+1):(3*(j-1)+3)] .= inv(WOOFER_CONFIG.INERTIA)*mpc_config.r_hat*RotZ(ψ)
			else
				mpc_config.B_c[7:12, (3*(j-1)+1):(3*(j-1)+3)] .= zeros(6, 3)
			end
		end

		# get the discretized B matrix via ZOH approximation:
		mpc_config.B_d .= mpc_config.dt*mpc_config.A_d*mpc_config.B_c
		mpc_config.B_lineq[12(i-1)+1:12(i-1)+12, 12(i-1)+1:12(i-1)+12] .= mpc_config.B_d
	end

	mpc_config.lb[1:12] .= -x0
	mpc_config.ub[1:12] .= -x0

	mpc_config.C_dense[13:12*(mpc_config.N+1), (12*mpc_config.N+13):(24*mpc_config.N+12)] .= mpc_config.B_lineq

	# update z_ref with state reference
	mpc_config.z_ref[1:12] .= x0
	mpc_config.z_ref[13:12*(mpc_config.N+1)] .= vec(x_ref)
	mpc_config.z_ref[12*(mpc_config.N+1)+1:24*(mpc_config.N)+12] .= repeat(mpc_config.u_ref, mpc_config.N)

	# add LQR cost to go to H
	mpc_config.H[12*(mpc_config.N)+1:12*(mpc_config.N)+12, 12*(mpc_config.N)+1:12*(mpc_config.N)+12] .= dare(	mpc_config.A_d, mpc_config.B_d,
																												mpc_config.Q, mpc_config.R)

	mpc_config.P .= sparse(2*mpc_config.H)
	mpc_config.q .= sparse(-2*mpc_config.H*mpc_config.z_ref)

	OSQP.setup!(mpc_config.prob; P=mpc_config.P, q=mpc_config.q, A=sparse(mpc_config.C_dense), l=mpc_config.lb, u=mpc_config.ub, verbose=false)
	results = OSQP.solve!(mpc_config.prob)

	mpc_config.z_result .= results.x

	forces .= (results.x)[12*(mpc_config.N+1)+1:12*(mpc_config.N+1)+12]

	println(forces)

end
