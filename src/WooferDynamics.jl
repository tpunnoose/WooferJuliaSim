function legJacobian(αᵢ::Vector{Float64}, abduction_offset = 0)
	# ported from WooferDynamics.py
	i = [1, 0, 0]
	j = [0, 1, 0]
	k = [0, 0, 1]

	beta = αᵢ[1]
	theta = αᵢ[2]
	r = αᵢ[3]

	# leg pointing straight down
	unrotated = [0, abduction_offset, -WOOFER_CONFIG.LEG_L + r]
	# vector from leg hub to foot in body coordinates
	p = RotYX(theta, beta)*unrotated

	# directions of positive joint movement
	theta_axis = RotX(beta)*j
	beta_axis = i
	radial_axis = k

	J = zeros(3,3)

	# dpdbeta
	J[:,1] .= cross(beta_axis, p)
	# dpdtheta
	J[:,2] .= cross(theta_axis, p)
	# dpdr
	J[:,3] .= RotYX(theta, beta)*radial_axis

	return J

end

function legJacobian2(α::Vector{Float64}, abduction_offset = 0)

	beta = α[1]
	theta = α[2]
	r = α[3]

	J = zeros(3,3)
	J[1,1] = 0
	J[1,2] = -(WOOFER_CONFIG.LEG_L + WOOFER_CONFIG.FOOT_RADIUS + r)*cos(theta)
	J[1,3] = -sin(theta)
	J[2,1] = (WOOFER_CONFIG.LEG_L + WOOFER_CONFIG.FOOT_RADIUS + r)*cos(beta)*cos(theta)
	J[2,2] = -(WOOFER_CONFIG.LEG_L + WOOFER_CONFIG.FOOT_RADIUS + r)*cos(beta)*sin(theta)
	J[2,3] = sin(beta)*cos(theta)
	J[3,1] = -(WOOFER_CONFIG.LEG_L + WOOFER_CONFIG.FOOT_RADIUS + r)*sin(beta)*cos(theta)
	J[3,2] = -(WOOFER_CONFIG.LEG_L + WOOFER_CONFIG.FOOT_RADIUS + r)*cos(beta)*sin(theta)
	J[3,3] = cos(beta)*cos(theta)

	return J
end

function force2Torque!(τ::Vector{Float64}, f::Vector{Float64}, α::Vector{Float64})
	τ[1:3] .= transpose(legJacobian(α[1:3]))*f[1:3]
	τ[4:6] .= transpose(legJacobian(α[4:6]))*f[4:6]
	τ[7:9] .= transpose(legJacobian(α[7:9]))*f[7:9]
	τ[10:12] .= transpose(legJacobian(α[10:12]))*f[10:12]
end

function forwardKinematics!(r_body::Vector{Float64}, α::Vector{Float64}, i::Int, abduction_offset=0)
	beta = α[1]
	theta = α[2]
	r = α[3]

	# leg pointing straight down
	unrotated = [0, abduction_offset, -WOOFER_CONFIG.LEG_L + r]
	# vector from leg hub to foot in body coordinates
	p = RotYX(theta, beta)*unrotated

	if i==1
		r_body .= p + [WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
	elseif i==2
		r_body .= p + [WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
	elseif i==3
		r_body .= p + [-WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
	else
		r_body .= p + [-WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
	end
end

function inverseKinematics!(α::Vector{T}, r_body::Vector{T}, i::Int, n::Int = 5, abduction_offset::Number=0) where {T<:Number}
	#=
	calculate joint angles via Newton's method
	α: initial guess (that will be updated in place)
	n: max number of iterations

	# TODO:
	=#

	r_guess = zeros(3)
	eps = 1e-2

	for k=1:n
		forwardKinematics!(r_guess, zeros(3), i)
		if norm(r_body - r_guess) < eps
			break
		end
		α .= pinv(legJacobian(α))* (r_body - r_guess) + α
	end

	print(norm(r_body - r_guess))

end
