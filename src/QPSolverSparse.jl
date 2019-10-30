"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

using OSQP
using ControlSystems
using BlockArrays
# using Parametron # TODO: use this

@with_kw struct MPCControllerParams
	# inverse of inertia tensor in body frame
	J_inv::Diagonal{Float64}

	# discretization length
	dt::Float64

	# planning horizon length
	N::Int64

	n::Int64 = 12
	m::Int64 = 12
	c::Int64 = 20

	d::Vector{Float64} = [0, 0, 0, 0, 0, 0, 0, 0, -9.81, 0, 0, 0]*dt

	Q::Diagonal{Float64}
	R::Diagonal{Float64}

	lb_friction::Vector{Float64}
	ub_friction::Vector{Float64}

	C_friction::Matrix{Float64}
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

	mpcConfig = MPCControllerParams(J_inv=J_inv, dt=dt, N=N, lb_friction=lb_friction, ub_friction=ub_friction, C_friction=C_i, Q=Q, R=R)
	return mpcConfig
end

function generateReferenceTrajectory!(x_ref::Array{T, 2}, x_curr::Vector{T}, x_des::Vector{T}, mpc_config::MPCControllerParams) where {T<:Number}
	# TODO: integrate the x,y,ψ position from the reference

	x_diff = 1/mpc_config.N * (x_des - x_curr)
	x_ref[:,1] .= x_curr
	for i in 1:mpc_config.N
		if i<mpc_config.N
			x_ref[:, i+1] .= x_curr + i*x_diff
		else
			x_ref[:, i+1] .= x_des
		end
	end
end

function solveFootForces!(forces::Vector{T}, x_ref::Array{T, 2}, contacts::Array{Int,2}, foot_locs::Array{T,2}, mpc_config::MPCControllerParams) where {T<:Number}
	# x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
	# contacts: 4xN+1 matrix of foot contacts over the planning horizon
	# foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon

	row_array = repeat([mpc_config.n, mpc_config.c], mpc_config.N)
	col_array = repeat([mpc_config.m, mpc_config.n], mpc_config.N)

	H = BlockArray(zeros(sum(col_array), sum(col_array)), col_array, col_array)
	C = BlockArray(zeros(sum(row_array), sum(col_array)), row_array, col_array)
	l = BlockArray(zeros(sum(row_array), 1), row_array, [1])
	u = BlockArray(zeros(sum(row_array), 1), row_array, [1])
	z_ref = BlockArray(zeros(sum(col_array), 1), col_array, [1])

	d = [0, 0, 0, 0, 0, 0, 0, 0, -9.81, 0, 0, 0]*dt
	u_ref = [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4
	# u_ref = zeros(12)

	A_c = zeros(12,12)
	B_c = zeros(12,12)
	A_d_i = zeros(12,12)
	B_d_i = zeros(12,12)
	V = zeros(12,12)

	G = zeros(3,3)
	cont_sys = zeros(24,24)
	r_hat = zeros(3,3)

	for i in 1:mpc_config.N
		## construct continuous time A, B matrix
		ψ_i = x_ref[6,i]

		for j in 1:4
			k = 3*(j-1)+1 # 1 -> 4 -> 7 -> 10
			if contacts[j,i] == 1
				B_c[7:9, k:k+2] .= 1/WOOFER_CONFIG.MASS*Matrix{Float64}(I, 3, 3)

				skewSymmetricMatrix!(r_hat, foot_locs[k:k+2, i])
				B_c[10:12, k:k+2] .= mpc_config.J_inv*r_hat*RotZ(ψ_i)
			else
				B_c[7:12, k:k+2] .= zeros(6, 3)
			end

			A_c[1:3, 7:9] = Matrix{Float64}(I, 3, 3)

			G[1,1] = sqrt(1 - 0.25*ψ_i^2)
			G[1,2] = 0.5*ψ_i
			G[2,1] = -0.5*ψ_i
			G[2,2] = sqrt(1 - 0.25*ψ_i^2)
			G[3,3] = sqrt(1 - 0.25*ψ_i^2)
			A_c[4:6, 10:12] .= G
		end

		cont_sys[1:12, 1:12] .= A_c
		cont_sys[1:12, 13:24] .= B_c
		disc_sys = exp(cont_sys*mpc_config.dt)

		# get the discretized A,B matrix via ZOH approximation:
		A_d_i = disc_sys[1:12, 1:12]
		B_d_i = disc_sys[1:12, 13:24]

		# put dynamics matrices into constraint matrix
		if i==1
			# set dynamics upper and lower bounds
			lu = -A_d_i*x_ref[:,1] - mpc_config.d
			lu = lu[:, :]
			setblock!(l, lu, i, 1)
			setblock!(u, lu, i, 1)
		else
			setblock!(C, A_d_i, 2*i-1, 2*(i-1))
			setblock!(l, -mpc_config.d[:,:], 2*i-1, 1)
			setblock!(u, -mpc_config.d[:,:], 2*i-1, 1)
		end
		setblock!(C, B_d_i, 2*i-1, 2*i-1)
		setblock!(C, -Matrix{Float64}(I, mpc_config.n, mpc_config.n), 2*i-1, 2*i)
		setblock!(C, mpc_config.C_friction, 2*i, 2*i-1)

		# set control friction and force bounds
		# TODO: only do this on first run?
		setblock!(l, mpc_config.lb_friction[:,:], 2*i, 1)
		setblock!(u, mpc_config.ub_friction[:,:], 2*i, 1)

		# populate reference trajectory
		setblock!(z_ref, u_ref[:,:], 2*i-1, 1)
		setblock!(z_ref, (x_ref[:, i+1])[:,:], 2*i, 1)

		# populate cost matrix
		setblock!(H, mpc_config.R, 2*i-1, 2*i-1)
		# TODO: test LQR cost to go if i==N
		if i==mpc_config.N
			V .= dare(A_d_i, B_d_i, Diagonal(mpc_config.Q), Diagonal(mpc_config.R))
			setblock!(H, V, 2*i, 2*i)
		else
			setblock!(H, mpc_config.Q, 2*i, 2*i)
		end
	end

	P = Array(2*H)
	q = Array(-2*H*z_ref)[:,1]
	l_qp = Array(l)[:,1]
	u_qp = Array(u)[:,1]

	prob = OSQP.Model()
	# something is screwy with the constraints
	OSQP.setup!(prob; P=sparse(P), q=q, A=sparse(C[1:12, :]), l=l_qp[1:12], u=u_qp[1:12], verbose=false)
	results = OSQP.solve!(prob)

	z_result = results.x

	forces .= (results.x)[1:12]

	println(forces)

end
