"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

using OSQP
using ControlSystems
using BlockArrays
using LinearAlgebra

using Convex
# using ECOS
# using SCS
# using Parametron # TODO: use this

include("WooferConfig.jl")

@with_kw struct MPCControllerParams
	# inverse of inertia tensor in body frame
	J_inv::Diagonal{Float64}

	# discretization length
	dt::Float64

	# planning horizon length
	N::Int64

	n::Int64 = 12
	m::Int64 = 12

	# Q::Diagonal{Float64} = Diagonal([1e3, 1e3, 1e2, 1e2, 1e2, 1e2, 1e3, 1e3, 1e3, 1, 1, 1e2])
	# R::Diagonal{Float64} = Diagonal([1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4])
	Q::Diagonal{Float64} = Diagonal([1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2])
	R::Diagonal{Float64} = Diagonal([1e-1, 1e-1, 1e-3, 1e-1, 1e-1, 1e-3, 1e-1, 1e-1, 1e-3, 1e-1, 1e-1, 1e-3])

	lb_friction::Vector{Float64}
	ub_friction::Vector{Float64}

	min_vert_force::Float64
	max_vert_force::Float64

	C_friction::Matrix{Float64}
	c::Int64 = size(C_friction)[1]

	row_array::Vector{Int64} = repeat([n, c], N)
	col_array::Vector{Int64} = repeat([m, n], N)

	# gravity affine term -> maps to vdot
	d::Vector{Float64} = [0, 0, 0, 0, 0, 0, 0, 0, -9.81, 0, 0, 0]*dt

	# block arrays
	H::BlockArray{Float64} = BlockArray(zeros(sum(col_array), sum(col_array)), col_array, col_array)
	C::BlockArray{Float64} = BlockArray(zeros(sum(row_array), sum(col_array)), row_array, col_array)
	l::BlockArray{Float64} = BlockArray(zeros(sum(row_array), 1), row_array, [1])
	u::BlockArray{Float64} = BlockArray(zeros(sum(row_array), 1), row_array, [1])
	z_ref::BlockArray{Float64} = BlockArray(zeros(sum(col_array), 1), col_array, [1])

	A_c::Matrix{Float64} = zeros(n,n)
	B_c::Matrix{Float64} = zeros(n,m)
	A_d_i::Matrix{Float64} = zeros(n,n)
	B_d_i::Matrix{Float64} = zeros(n,m)
	V::Matrix{Float64} = zeros(n,n)

	G::Matrix{Float64} = zeros(3,3)
	cont_sys::Matrix{Float64} = zeros(n+m,n+m)
	disc_sys::Matrix{Float64} = zeros(n+m,n+m)
	r_hat::Matrix{Float64} = zeros(3,3)

	z_result = zeros(sum(col_array))

	prob::OSQP.Model
end

function initMPCControllerConfig(dt::Number, N::Integer)
	# z = (x0, ..., xN, u0, ..., u0, ..., uN)
	J_inv = inv(WOOFER_CONFIG.INERTIA)

	prob = OSQP.Model()

	mu = 1.5
	min_vert_force = 1
	max_vert_force = 133

	C_friction = zeros(20, 12)

	# friction constraint matrix for each u_i
	for i in 1:4
		i_3 = 3*(i-1)
		i_5 = 5*(i-1)
		# fz >= 0
		C_friction[i_5+1, i_3+3] = 1
		# ufz+fx >= 0
		C_friction[i_5+2, i_3+1] = 1
		C_friction[i_5+2, i_3+3] = mu
		# fx-ufz <= 0
		C_friction[i_5+3, i_3+1] = 1
		C_friction[i_5+3, i_3+3] = -mu
		# ufz+fy >= 0
		C_friction[i_5+4, i_3+2] = 1
		C_friction[i_5+4, i_3+3] = mu
		# fy-ufz <= 0
		C_friction[i_5+5, i_3+2] = 1
		C_friction[i_5+5, i_3+3] = -mu
	end

	# C_friction = Matrix{Float64}(I, 20, 12)

	lb_friction = zeros(20)
	ub_friction = zeros(20)

	# friction lower and upper bounds for each u_i
	for i in 1:5:20
		# fz >= 0
		lb_friction[i] = min_vert_force
		ub_friction[i] = max_vert_force
		# ufz+fx >= 0
		lb_friction[i+1] = 0
		ub_friction[i+1] = 1e5
		# fx-ufz <= 0
		lb_friction[i+2] = -1e5
		ub_friction[i+2] = 0
		# ufz+fy >= 0
		lb_friction[i+3] = 0
		ub_friction[i+3] = 1e5
		# fy-ufz >= 0
		lb_friction[i+4] = -1e5
		ub_friction[i+4] = 0
	end

	Q = Diagonal([1e3, 1e3, 5e4, 1e2, 1e2, 1e2, 1e3, 1e3, 1e3, 1, 1, 1e2])
	R = Diagonal([1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4])

	mpcConfig = MPCControllerParams(J_inv=J_inv, dt=dt, N=N, lb_friction=lb_friction,
									ub_friction=ub_friction, C_friction=C_friction, prob=prob,
									min_vert_force=min_vert_force, max_vert_force=max_vert_force)
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

function solveFootForces!(forces::Vector{T}, x_ref::Array{T, 2}, contacts::Array{Int,2}, foot_locs::Array{T,2}, mpc_config::MPCControllerParams, use_lqr::Bool=false) where {T<:Number}
	# x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
	# contacts: 4xN+1 matrix of foot contacts over the planning horizon
	# foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon

	if use_lqr
		u_ref = [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4
	else
		u_ref = zeros(12)
	end

	for i in 1:mpc_config.N
		## construct continuous time A, B matrix
		ψ_i = x_ref[6,i]

		for j in 1:4
			k = 3*(j-1)+1 # 1 -> 4 -> 7 -> 10
			if contacts[j,i] == 1
				mpc_config.B_c[7:9, k:k+2] .= 1/WOOFER_CONFIG.MASS*Matrix{Float64}(I, 3, 3)

				skewSymmetricMatrix!(mpc_config.r_hat, foot_locs[k:k+2, i])
				mpc_config.B_c[10:12, k:k+2] .= mpc_config.J_inv*mpc_config.r_hat*RotZ(ψ_i)
				for i=1:5:20
					mpc_config.lb_friction[i] = mpc_config.min_vert_force
					mpc_config.ub_friction[i] = mpc_config.max_vert_force
				end
			else
				mpc_config.B_c[7:12, k:k+2] .= zeros(6, 3)
				for i=1:5:20
					mpc_config.lb_friction[i] = 0
					mpc_config.ub_friction[i] = 0
				end
			end

			mpc_config.A_c[1:3, 7:9] = Matrix{Float64}(I, 3, 3)

			mpc_config.G[1,1] = sqrt(1 - 0.25*ψ_i^2)
			mpc_config.G[1,2] = 0.5*ψ_i
			mpc_config.G[2,1] = -0.5*ψ_i
			mpc_config.G[2,2] = sqrt(1 - 0.25*ψ_i^2)
			mpc_config.G[3,3] = sqrt(1 - 0.25*ψ_i^2)
			mpc_config.A_c[4:6, 10:12] .= mpc_config.G
		end

		mpc_config.cont_sys[1:12, 1:12] .= mpc_config.A_c
		mpc_config.cont_sys[1:12, 13:24] .= mpc_config.B_c
		mpc_config.disc_sys .= exp(mpc_config.cont_sys*mpc_config.dt)

		# get the discretized A,B matrix via ZOH approximation:
		mpc_config.A_d_i .= mpc_config.disc_sys[1:12, 1:12]
		mpc_config.B_d_i .= mpc_config.disc_sys[1:12, 13:24]

		# put dynamics matrices into constraint matrix
		if i==1
			# set dynamics upper and lower bounds
			lu = -mpc_config.A_d_i*x_ref[:,1] - mpc_config.d
			lu = lu[:, :]
			setblock!(mpc_config.l, lu, i, 1)
			setblock!(mpc_config.u, lu, i, 1)
		else
			setblock!(mpc_config.C, mpc_config.A_d_i, 2*i-1, 2*(i-1))
			setblock!(mpc_config.l, -mpc_config.d[:,:], 2*i-1, 1)
			setblock!(mpc_config.u, -mpc_config.d[:,:], 2*i-1, 1)
		end
		setblock!(mpc_config.C, mpc_config.B_d_i, 2*i-1, 2*i-1)
		setblock!(mpc_config.C, -Matrix{Float64}(I, mpc_config.n, mpc_config.n), 2*i-1, 2*i)
		setblock!(mpc_config.C, mpc_config.C_friction, 2*i, 2*i-1)

		# set control friction and force bounds
		# TODO: only do this on first run?
		setblock!(mpc_config.l, mpc_config.lb_friction[:,:], 2*i, 1)
		setblock!(mpc_config.u, mpc_config.ub_friction[:,:], 2*i, 1)

		# populate reference trajectory
		setblock!(mpc_config.z_ref, u_ref[:,:], 2*i-1, 1)
		setblock!(mpc_config.z_ref, (x_ref[:, i+1])[:,:], 2*i, 1)

		# populate cost matrix
		setblock!(mpc_config.H, mpc_config.R, 2*i-1, 2*i-1)
		# TODO: test LQR cost to go if i==N
		if use_lqr && (i==mpc_config.N)
			mpc_config.V .= dare(mpc_config.A_d_i, mpc_config.B_d_i, Diagonal(mpc_config.Q), Diagonal(mpc_config.R))
			# symmetrize V -> shouldn't have to do this
			mpc_config.V .= 0.5*mpc_config.V + 0.5*mpc_config.V'
			setblock!(mpc_config.H, mpc_config.V, 2*i, 2*i)
		else
			setblock!(mpc_config.H, mpc_config.Q, 2*i, 2*i)
		end
	end

	P = Array(2*mpc_config.H)
	q = Array(-2*mpc_config.H*mpc_config.z_ref)[:,1]
	l_qp = Array(mpc_config.l)[:,1]
	u_qp = Array(mpc_config.u)[:,1]

	prob = OSQP.Model()
	OSQP.setup!(prob; P=sparse(P), q=q, A=sparse(mpc_config.C[1:12,:]), l=l_qp[1:12], u=u_qp[1:12], verbose=false, polish=1, eps_abs=1e-12, eps_rel=1e-12)
	results = OSQP.solve!(prob)
	mpc_config.z_result .= results.x

	# z = Variable(mpc_config.N*(mpc_config.m + mpc_config.n))
	# problem = minimize(quadform(z - mpc_config.z_ref, P), u_qp >= mpc_config.C*z, mpc_config.C*z >= l_qp)
	# solve!(problem, ECOSSolver(maxit=1000))
	# mpc_config.z_result .= evaluate(z)[:,1]

	forces .= (mpc_config.z_result)[1:12]


end
