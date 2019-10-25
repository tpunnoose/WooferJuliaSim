using ControlSystems
using OSQP

@with_kw struct QPParams
	V::Array{Float64, 2}
	u0::Vector{Float64}
	x0::Vector{Float64}
	dt::Float64
	A_d::Array{Float64, 2}
	B_d::Array{Float64, 2}
	Q::Array{Float64, 2}
	R::Array{Float64, 2}
	eps::Float64 = 0.75
end

function initQPParams(dt::AbstractFloat, x0::Vector, ψ::AbstractFloat = 0.0)
	r1 = zeros(3)
    r2 = zeros(3)
    r3 = zeros(3)
    r4 = zeros(3)

	forwardKinematics!(r1, zeros(3), 1)
	forwardKinematics!(r2, zeros(3), 2)
	forwardKinematics!(r3, zeros(3), 3)
	forwardKinematics!(r4, zeros(3), 4)

	foot_locs = zeros(12)

	foot_locs[1:3] .= r1
	foot_locs[4:6] .= r2
	foot_locs[7:9] .= r3
	foot_locs[10:12] .= r4

	A_c = zeros(12,12)
	A_c[1:3, 7:9] = Matrix{Float64}(I, 3, 3)

	G = zeros(3,3)
	G[1,1] = sqrt(1 - 0.25*ψ^2)
	G[1,2] = 0.5*ψ
	G[2,1] = -0.5*ψ
	G[2,2] = sqrt(1 - 0.25*ψ^2)
	G[3,3] = sqrt(1 - 0.25*ψ^2)
	A_c[4:6, 10:12] .= G

	B_c = zeros(12,12)
	r_hat = zeros(3,3)

	for j in 1:4
		B_c[7:9, (3*(j-1)+1):(3*(j-1)+3)] .= 1/WOOFER_CONFIG.MASS*Matrix{Float64}(I, 3, 3)

		skewSymmetricMatrix!(r_hat, foot_locs[(3*(j-1)+1):(3*(j-1)+3)])
		B_c[10:12, (3*(j-1)+1):(3*(j-1)+3)] .= inv(WOOFER_CONFIG.INERTIA)*r_hat*RotZ(ψ)
	end

	# matrix exponential approximation
	# A_d = (Matrix{Float64}(I, 10, 10) + A_c*dt + (A_c*dt)^2/2 + (A_c*dt)^3/6 + (A_c*dt)^4/24 + (A_c*dt)^5/120 + (A_c*dt)^6/720)
	# B_d = dt*(Matrix{Float64}(I, 10, 10) + A*dt + (A*dt)^2/2 + (A*dt)^3/6 + (A*dt)^4/24 + (A*dt)^5/120 + (A*dt)^6/720)*B_c

	cont_sys = zeros(24,24)
	cont_sys[1:12, 1:12] .= A_c
	cont_sys[1:12, 13:24] .= B_c
	disc_sys = exp(cont_sys*dt)

	A_d = disc_sys[1:12, 1:12]
	B_d = disc_sys[1:12, 13:24]

	# # calculate u around new point
	# x1 = [0.34, 0, 0, 0.1, -0.1, 0, 0, 0, 0]
	# println(B_d\((Matrix{Float64}(I, 9, 9) - A_d)*x1 + dt*[0.0, 0, 0, 0.0, -0.0, -9.81, 0, 0, 0]))

	# Q = Diagonal([1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e3, 1e3, 1e3, 1, 1, 1e2])
	Q = Diagonal([1e3, 1e3, 1e2, 1e2, 1e2, 1e2, 1e3, 1e3, 1e3, 1, 1, 1e2])
	R = Diagonal([1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4])

	#
	# println(rank(ctrb(A_d, B_d)))

	V = dare(A_d, B_d, Q, R)
	# x0 = [0.34, 0, 0, 0, 0, 0, 0, 0, 0]
	# u0 = [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4

	# calculate the steady state control for the desired state
	u0 = pinv(B_d)*((Matrix{Float64}(I, 12, 12) - A_d)*x0 + dt*[0.0, 0, 0, 0.0, 0.0, 0, 0, 0, 9.81, 0, 0, 0])

	return QPParams(x0=x0, u0=u0, dt=dt, A_d=A_d, B_d=B_d, Q=Q, R=R, V=V)
end

function qpBalance!(forces::Vector{T}, x::Vector{T}, params::QPParams) where {T<:Number}
	# LQR solved via (dense) one step QP:
	P = 2*(params.R + params.B_d'*params.V*params.B_d)
	q = (2*(x - params.x0)'*params.A_d'*params.V*params.B_d)'


	# TODO: put C/lb/ub in the QPParams structure
	C = zeros(20, 12)
	mu = 0.6
	min_vert_force = 1
	max_vert_force = 133

	for i in 1:4
		# fz >= 0
		C[(i-1)*5+1, (i-1)*3+3] = 1
		# ufz+fx >= 0
		C[(i-1)*5+2, (i-1)*3+1] = 1
		C[(i-1)*5+2, (i-1)*3+3] = mu
		# fx-ufz <= 0
		C[(i-1)*5+3, (i-1)*3+1] = 1
		C[(i-1)*5+3, (i-1)*3+3] = -mu
		# ufz+fy >= 0
		C[(i-1)*5+4, (i-1)*3+2] = 1
		C[(i-1)*5+4, (i-1)*3+3] = mu
		# fy-ufz <= 0
		C[(i-1)*5+5, (i-1)*3+2] = 1
		C[(i-1)*5+5, (i-1)*3+3] = -mu
	end

	lb = zeros(20)
	ub = zeros(20)

	for i in 1:4
		# fz >= 0
		lb[(i-1)*5+1] = min_vert_force
		ub[(i-1)*5+1] = max_vert_force
		# ufz+fx >= 0
		lb[(i-1)*5+2] = 0
		ub[(i-1)*5+2] = Inf
		# fx-ufz <= 0
		lb[(i-1)*5+3] = -Inf
		ub[(i-1)*5+3] = 0
		# ufz+fy >= 0
		lb[(i-1)*5+4] = 0
		ub[(i-1)*5+4] = Inf
		# fy-ufz >= 0
		lb[(i-1)*5+5] = -Inf
		ub[(i-1)*5+5] = 0
	end

	lb = lb - C*params.u0
	ub = ub - C*params.u0

	# TODO: add friction cone constraints, add warm starting
	prob = OSQP.Model()
	OSQP.setup!(prob; P=sparse(P), q=q, A=sparse(C), l=lb, u=ub, verbose=false)
	# OSQP.setup!(prob; P=sparse(P), q=q, verbose=false)
	results = OSQP.solve!(prob)

	forces .= results.x + params.u0
end

function skewSymmetricMatrix!(A::Matrix, a::Vector)
	A[1,2] = -a[3]
	A[1,3] = a[2]
	A[2,1] = a[3]
	A[2,3] = -a[1]
	A[3,1] = -a[2]
	A[3,2] = a[1]
end
