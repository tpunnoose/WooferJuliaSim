using ControlSystems

@with_kw struct LQRParams
	L::Array{Float64, 2}
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

function initLQRParams(dt::AbstractFloat, x0::Vector, ψ::AbstractFloat = 0.0)
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
		println(foot_locs[(3*(j-1)+1):(3*(j-1)+3)])
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

	Q = Diagonal([1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e3, 1e3, 1e3, 1, 1, 1e2])
	R = Diagonal([1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2, 1e-4])
	#
	# println(A_d)
	# println(B_d)
	#
	# println(rank(ctrb(A_d, B_d)))

	L = dlqr(A_d, B_d, Q, R)
	V = dare(A_d, B_d, Q, R)
	# x0 = [0.34, 0, 0, 0, 0, 0, 0, 0, 0]
	u0 = [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4

	return LQRParams(L=L, x0=x0, u0=u0, dt=dt, A_d=A_d, B_d=B_d, Q=Q, R=R, V=V)
end

function lqrBalance!(forces::Vector{T}, x::Vector{T}, joint_pos::Vector{T}, lqr_params::LQRParams) where {T<:Number}
	forces .= lqr_params.u0 - lqr_params.L*(x - lqr_params.x0)
	println(forces)
end

function skewSymmetricMatrix!(A::Matrix, a::Vector)
	A[1,2] = -a[3]
	A[1,3] = a[2]
	A[2,1] = a[3]
	A[2,3] = -a[1]
	A[3,1] = -a[2]
	A[3,2] = a[1]
end
