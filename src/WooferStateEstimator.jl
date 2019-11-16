@with_kw mutable struct StateEstimatorParams
	# current estimate
	x::Vector{Float64}
	P::Matrix{Float64}

	xdot::Vector{Float64} = zeros(size(x)[1])
	x_::Vector{Float64} = zeros(size(x)[1])
	x_plus::Vector{Float64} = zeros(size(x)[1])
	P_::Matrix{Float64} = zeros(size(x)[1], size(x)[1])
	P_plus::Matrix{Float64} = zeros(size(x)[1], size(x)[1])

	A::Matrix{Float64} = zeros(size(x)[1], size(x)[1])
	C::Matrix{Float64} = zeros(4, size(x)[1])
	S::Matrix{Float64} = zeros(4,4)
	y_meas::Vector{Float64} = zeros(4)
	y_pred::Vector{Float64} = zeros(4)
	nu::Vector{Float64} = zeros(4)
	K::Matrix{Float64} = zeros(size(x)[1],4)
	Q::Matrix{Float64}
	R::Matrix{Float64}
	g_n::Vector{Float64} = [0, 0, 9.81]
	g_b::Vector{Float64} = [0, 0, 9.81]
	qs::Float64 = 0.0
	qv::Vector{Float64} = [0, 0, 0]
	dqvdphi::Vector{Float64} = [0.5, 0, 0]
	dqvdtheta::Vector{Float64} = [0, 0.5, 0]
	dqsdphi::Float64 = 0.0
	dqsdtheta::Float64 = 0.0

	r_rel::Vector{Float64} = zeros(3)

	dt::Float64

	r_fr::Vector{Float64} = [WOOFER_CONFIG.LEG_FB, 	-WOOFER_CONFIG.LEG_LR, 0]
	r_fl::Vector{Float64} = [WOOFER_CONFIG.LEG_FB, 	 WOOFER_CONFIG.LEG_LR, 0]
	r_br::Vector{Float64} = [-WOOFER_CONFIG.LEG_FB, 	-WOOFER_CONFIG.LEG_LR, 0]
	r_bl::Vector{Float64} = [-WOOFER_CONFIG.LEG_FB, 	 WOOFER_CONFIG.LEG_LR, 0]
end

function skewSymmetricMatrix(q::Vector{Float64})
	return [0.0 -q[3] q[2]; q[3] 0.0 -q[1]; -q[2] q[1] 0.0]
end

function L_q(q::Vector{Float64})
	return [q[1] -q[2:4]';q[2:4] q[1]*I+skewSymmetricMatrix(q[2:4])]
end

function R_q(q::Vector{Float64})
	return [q[1] -q[2:4]';q[2:4] q[1]*I-skewSymmetricMatrix(q[2:4])]
end

function initStateEstimator(dt::T, x::Vector{T}, P::Diagonal{T}, Q::Diagonal{T}, R::Diagonal{T}) where {T<:Number} # default Q, R?
	return StateEstimatorParams(dt=dt, x=x, P=P, Q=Q, R=R)
end

function stateEstimatorUpdate(dt::AbstractFloat, r_hat::Vector, q_hat::Vector, joint_pos::Vector, joint_vel::Vector, est_params::StateEstimatorParams, torques::Vector)
	x_r = est_params.x[1:3]
	x_ϕ = est_params.x[4:6]
	x_v = est_params.x[7:9]
	x_w = est_params.x[10:12]

	m = WOOFER_CONFIG.MASS
	τ_F_const = 0.005

	# placeholder
	J = Matrix{Float64}(I,3,3)

	V = zeros(3,4)
	V[:,2:4] = Matrix{Float64}(I,3,3)

	A_t = zeros(12,12)
	A_t[1:3,7:9] = Matrix{Float64}(I, 3, 3)
	A_t[4:6,4:6] = 0.5*V*R_q([0;x_w])*[x_ϕ'/sqrt(abs(1-x_ϕ'*x_ϕ));I]
	A_t[4:6,10:12] = 0.5*(L_q([1;x_ϕ])*V')[2:4,:]
	A_t[10:12,10:12] = inv(J)*(skewSymmetricMatrix(x_w)*J)

	B_t = zeros(12,6)
	B_t[7:9,1:3] = Matrix{Float64}(I,3,3)/m
	B_t[10:12,4:6] = inv(J)

	C_t = [I zeros(3,9); zeros(3,3) (V*L_q([1;q_hat])')[:,2:4] zeros(3,6)]

	y_hat = C_t * est_params.x[1:12]
	y_true = [r_hat; q_hat]

	F = zeros(3)
	τ = zeros(3)
	leg_vec = zeros(3)

	for i in 1:4
		# forwardKinematics!(est_params.r_rel, joint_pos[3*(i-1)+1:3*(i-1)+3], i)
		# if i==1
		# 	leg_vec = est_params.r_rel - [WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
		# elseif i==2
		# 	leg_vec = est_params.r_rel - [WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
		# elseif i==3
		# 	leg_vec = est_params.r_rel - [-WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
		# else
		# 	leg_vec = est_params.r_rel - [-WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
		# end

		forwardKinematics!(leg_vec, joint_pos[3*(i-1)+1:3*(i-1)+3], i)

		leg_torque = torques[3*(i-1)+1:3*(i-1)+3]
		leg_hat = skewSymmetricMatrix(leg_vec)
		F_leg = zeros(3)
		F_leg[1:2] = τ_F_const*(pinv(leg_hat) * leg_torque)[1:2]
		if i==1 || i==2
			F_leg[3] = -τ_F_const*leg_torque[3]
		else
			F_leg[3] = τ_F_const*leg_torque[3]
		end
		τ_leg = τ_F_const*cross(est_params.r_rel, F_leg)

		F = F + F_leg

		if i==2 || i==4
			τ = [τ[1] + τ_leg[1], τ[2] - τ_leg[2], τ[3] + τ_leg[3]]
		else
			τ = τ + τ_leg
		end

	end

	u = [F; τ]
	# est_params.x_ = A_t*est_params.x[1:12]
	# forward prop dynamics
	est_params.x_[1:3] = est_params.x_[1:3] + dt * est_params.x_[7:9]
	est_params.x_[4:6] = (0.5 * L_q([1;x_ϕ]) * V' * x_w)[2:4]
	est_params.x_[7:9] = F/m
	est_params.x_[10:12] = inv(J) * (τ - cross(x_w, J*x_w))


	est_params.P_ = A_t*est_params.P*A_t' + est_params.Q
	K_s = est_params.P_ * C_t' * inv(C_t * est_params.P_ * C_t' + est_params.R)
	est_params.x_plus = est_params.x_ + K_s*(y_true - y_hat)
	est_params.P_plus = est_params.P_ - K_s*C_t*est_params.P_

	est_params.x = est_params.x_plus
	est_params.P = est_params.P_plus

end

function stateEstimatorUpdate2(dt::AbstractFloat, a_b::Vector, om_b::Vector, joint_pos::Vector, joint_vel::Vector, contacts::Vector, est_params::StateEstimatorParams)
	"""
	EKF State Estimator
	x = [z_b, phi, theta, v_b, b_a, b_om]

	TODO: need to make everything in place (skewSymmetricMatrix -> skewSymmetricMatrix!)
	"""
	@show est_params.x
	# scalar part of quaternion
	est_params.qs = sqrt(1 - 0.25*(est_params.x[2]^2 + est_params.x[3]^2))

	# vector part of quaternion
	est_params.qv[1] = 0.5*est_params.x[2]
	est_params.qv[2] = 0.5*est_params.x[3]

	# calculate gravity vector rotated into body frame (from current estimate)
	est_params.g_b = 	est_params.g_n + 2*skewSymmetricMatrix(est_params.qv) *
						(skewSymmetricMatrix(est_params.qv)*est_params.g_n - est_params.qs*est_params.g_n)

	# d/dt z_b
	est_params.xdot[1] = est_params.x[6]
	# d/dt phi
	est_params.xdot[2] = (om_b[1] - est_params.x[10]) * est_params.qs +
							0.5*est_params.x[3]*(om_b[3] - est_params.x[12])
	# d/dt theta
	est_params.xdot[3] = (om_b[2] - est_params.x[11]) * est_params.qs -
							0.5*est_params.x[2]*(om_b[3] - est_params.x[12])
	# d/dt v_b
	est_params.xdot[4:6] = (a_b - est_params.x[7:9]) - est_params.g_b -
							skewSymmetricMatrix(om_b - est_params.x[10:12])*est_params.x[4:6]
	# no predicted changes in biases
	est_params.xdot[7:12] .= 0

	est_params.dqsdphi = -0.25*est_params.x[2]/est_params.qs
	est_params.dqsdtheta = -0.25*est_params.x[3]/est_params.qs

	# ∂z_b/∂z_b
	est_params.A[1,1] = 1
	# ∂z_b/∂v_b
	est_params.A[1,6] = dt

	# ∂phi/∂phi
	est_params.A[2,2] = 1 - 0.25*(om_b[1] - est_params.x[10])*est_params.x[2]*
							dt/est_params.qs
	# ∂phi/∂theta
	est_params.A[2,3] = (-0.25*(om_b[1] - est_params.x[10])*est_params.x[3]/
							est_params.qs + 0.5*(om_b[3] - est_params.x[12]))*dt
	# ∂phi/∂b_w
	est_params.A[2,10] = -est_params.qs*dt
	est_params.A[2,12] = -0.5*est_params.x[3]*dt

	# ∂theta/∂phi
	est_params.A[3,2] = (-0.25*(om_b[2] - est_params.x[11])*est_params.x[2]/
							est_params.qs - 0.5*(om_b[3] - est_params.x[12]))*dt
	# ∂theta/∂theta
	est_params.A[3,3] = 1 - 0.25*(om_b[2] - est_params.x[11])*est_params.x[3]*
							dt/est_params.qs
	# ∂theta/∂b_w
	est_params.A[3,11] = -est_params.qs*dt
	est_params.A[3,12] = 0.5*est_params.x[2]*dt

	# ∂v_b/∂phi
	est_params.A[4:6, 2] = -dt*((-2*skewSymmetricMatrix(est_params.qv) * skewSymmetricMatrix(est_params.g_n) -
				 			2*skewSymmetricMatrix(skewSymmetricMatrix(est_params.qv) * est_params.g_n) +
							2*est_params.qs*skewSymmetricMatrix(est_params.g_n)) * est_params.dqvdphi -
							2*est_params.dqsdphi*skewSymmetricMatrix(est_params.qv) * est_params.g_n)
	# ∂v_b/∂theta
	est_params.A[4:6, 3] = -dt*((-2*skewSymmetricMatrix(est_params.qv) * skewSymmetricMatrix(est_params.g_n) -
				 			2*skewSymmetricMatrix(skewSymmetricMatrix(est_params.qv) * est_params.g_n) +
							2*est_params.qs*skewSymmetricMatrix(est_params.g_n)) * est_params.dqvdtheta -
							2*est_params.dqsdtheta*skewSymmetricMatrix(est_params.qv) * est_params.g_n)
	# ∂v_b/∂v_b
	est_params.A[4:6, 4:6] .= Matrix{Float64}(I, 3, 3) - skewSymmetricMatrix(om_b - est_params.x[10:12])*dt
	# ∂v_b/∂b_a
	est_params.A[4:6, 7:9] .= -Matrix{Float64}(I, 3, 3)*dt
	# ∂v_b/∂b_w
	est_params.A[4:6, 10:12] .= -skewSymmetricMatrix(est_params.x[4:6])*dt

	# ∂b_a/∂b_a
	est_params.A[7:9, 7:9] .= Matrix{Float64}(I, 3, 3)

	# ∂b_w/∂b_w
	est_params.A[10:12, 10:12] .= Matrix{Float64}(I, 3, 3)

	# prediction step

	est_params.x_ = est_params.x + est_params.xdot*dt
	est_params.P_ = est_params.A * est_params.P * est_params.A' + est_params.Q

	# sequential measurement update
	est_params.x_plus = est_params.x_
	est_params.P_plus = est_params.P_

	for i in 1:4
		# only update the feet that are in contact
		if contacts[i] == 1
			forwardKinematics!(est_params.r_rel, joint_pos[3*(i-1)+1:3*(i-1)+3], i)
			est_params.C[1,1] = 1
			est_params.C[2:4,4:6] = Matrix{Float64}(I, 3, 3)
			est_params.C[2:4,10:12] = skewSymmetricMatrix(est_params.r_rel)

			# assuming that the foot is on level ground
			est_params.y_meas[1] = 0 - est_params.r_rel[3]
			est_params.y_meas[2:4] = -legJacobian(joint_pos[3*(i-1)+1:3*(i-1)+3])*joint_vel[3*(i-1)+1:3*(i-1)+3]

			est_params.y_pred[1] = est_params.x_plus[1]
			est_params.y_pred[2:4] = (est_params.x_plus[4:6] + skewSymmetricMatrix(om_b - est_params.x_plus[10:12])*est_params.r_rel)

			est_params.nu = est_params.y_meas - est_params.y_pred
			est_params.S = est_params.C*est_params.P_plus*est_params.C' + est_params.R

			est_params.K = est_params.P_plus * est_params.C' * inv(est_params.S)
			# println((est_params.nu' * inv(est_params.S) * est_params.nu))

			# outlier detection here for slipping feet
			if (est_params.nu' * inv(est_params.S) * est_params.nu) < Inf
				est_params.x_plus = est_params.x_plus + est_params.K * est_params.nu
				est_params.P_plus = est_params.P_plus - est_params.K * est_params.S * est_params.K'
			end
		end
	end

	est_params.x = est_params.x_plus
	est_params.P = est_params.P_plus
end
