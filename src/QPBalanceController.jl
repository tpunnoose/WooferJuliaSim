using OSQP

@with_kw struct BalanceControllerConfig
	wn_cart::Float64 = 20
	zeta_cart::Float64 = 0.8
	kp_cart::Float64 = wn_cart^2
	kd_cart::Float64 = 2*wn_cart*zeta_cart

	wn_ang::Float64 = 20
	zeta_ang::Float64 = 0.8
	kp_ang::Float64 = wn_ang^2
	kd_ang::Float64 = 2*wn_ang*zeta_ang

	mu::Float64 = 0.75
	beta::Float64 = 5.0e-2
	gamma::Float64 = 2.0e2
	alpha::Float64 = 1.0e-3

	C::SparseMatrixCSC{Float64, Int64}
	lb::Vector{Float64}
	ub::Vector{Float64}

	active_feet::Vector{Int64}
end

function initBalanceController()
	C = zeros(20, 12)

	# FIXME: mu, active_feet is really defined here and not in the config struct
	mu = 1.0
	active_feet = [1, 1, 1, 1]
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

	C = sparse(C)

	lb = zeros(20)
	ub = zeros(20)

	for i in 1:4
		# fz >= 0
		lb[(i-1)*5+1] = min_vert_force*active_feet[i]
		ub[(i-1)*5+1] = max_vert_force*active_feet[i]
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

	config = BalanceControllerConfig(mu=mu, C=C, lb=lb, ub=ub, active_feet=active_feet)

	return config
end

function balanceController!(torques, x, joint_pos, p_ref, o_ref, f_prev, config::BalanceControllerConfig)
	# position PD controller
	a_d_com = config.kp_cart*(p_ref - x[1:3]) + config.kd_cart*(-x[7:9]) + [0, 0, 9.81]
	f_d	= WOOFER_CONFIG.MASS * a_d_com

	# orientation PD controller
	α_d_com = config.kp_ang*(o_ref - x[4:6]) + config.kd_ang*(-x[10:12])
	τ_d = WOOFER_CONFIG.INERTIA * α_d_com

	ref_wrench = zeros(6)
	ref_wrench[1:3] .= f_d
	ref_wrench[4:6] .= τ_d

	# active feet
	active_feet = [1, 1, 1, 1]

	foot_locs = zeros(12)
	r_body = zeros(3)

	R = zeros(3,3)
	CRP2DCM!(R, x[4:6])

	for i in 1:4
		if active_feet[i] == 1
			forwardKinematics!(r_body, joint_pos[3*(i-1)+1:3*(i-1)+3], i)
			foot_locs[3*(i-1)+1:3*(i-1)+3] .= R * r_body
		end
	end

	# print("Foot locations: ")
	# println(foot_locs)

	f_world = zeros(12)
	solveFootForces!(f_world, ref_wrench, active_feet, foot_locs, f_prev, config)

	for i in 1:4
		if active_feet[i] == 1
			f_body = -transpose(R) * f_world[3*(i-1)+1:3*(i-1)+3]
			torques[3*(i-1)+1:3*(i-1)+3] .= transpose(legJacobian(joint_pos[3*(i-1)+1:3*(i-1)+3])) * f_body
		end
	end

end

function skewSymmetricMatrix(a::Vector)
	a_hat = zeros(3,3)
	a_hat[1,2] = -a[3]
	a_hat[1,3] = a[2]
	a_hat[2,1] = a[3]
	a_hat[2,3] = -a[1]
	a_hat[3,1] = -a[2]
	a_hat[3,2] = a[1]

	return a_hat
end

function solveFootForces!(forces, ref_wrench, active_feet, foot_locs, f_prev, config::BalanceControllerConfig)
	A = zeros(6,12)
	A[1:3, 1:3] .= Matrix{Float64}(I, 3, 3)
	A[1:3, 4:6] .= Matrix{Float64}(I, 3, 3)
	A[1:3, 7:9] .= Matrix{Float64}(I, 3, 3)
	A[1:3, 10:12] .= Matrix{Float64}(I, 3, 3)

	A[4:6, 1:3] .= skewSymmetricMatrix(foot_locs[1:3])
	A[4:6, 4:6] .= skewSymmetricMatrix(foot_locs[4:6])
	A[4:6, 7:9] .= skewSymmetricMatrix(foot_locs[7:9])
	A[4:6, 10:12] .= skewSymmetricMatrix(foot_locs[10:12])

	K = Diagonal([1e-1, 1e-1, 1e-1, config.gamma, config.gamma, config.gamma])

	P = 2*(transpose(A)*K*A + (config.alpha + config.beta)*Matrix{Float64}(I, 12, 12))
	P = sparse(P)

	q = -2*(transpose(A)*K*ref_wrench + config.beta*f_prev)

	prob = OSQP.Model()

	# TODO: not setup every time, only update?
	OSQP.setup!(prob; P=P, q=q, A=config.C, l=config.lb, u=config.ub, verbose=false)
	results = OSQP.solve!(prob)

	forces .= results.x
end
