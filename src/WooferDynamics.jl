using ForwardDiff

function forwardKinematics(α::Vector{T}) where {T<:Real}
	beta = α[1]
	ϕ_1 = α[2]
	ϕ_2 = -α[3]


	# coordinates of leg pointing straight down
	γ = 0.5*(pi - ϕ_1 - ϕ_2)
	θ = pi/2 - γ - ϕ_1

	# length of leg:
	d = WOOFER_CONFIG.l0*sin(γ)
	h1 = WOOFER_CONFIG.l0*cos(γ)
	h2 = sqrt(WOOFER_CONFIG.l1^2 - d^2)
	L = h1 + h2

	unrotated = [L*sin(θ), WOOFER_CONFIG.ABDUCTION_OFFSET, -L*cos(θ)]

	# rotated by abduction
	return RotX(beta) * unrotated
end

function forwardKinematics!(r_body::Vector{T}, α::Vector{T}, i::Int64) where {T<:Real}
	beta = α[1]
	ϕ_1 = α[2] # TODO: check the sign of these
	ϕ_2 = α[3] # TODO: check the sign of these


	# coordinates of leg pointing straight down
	γ = 0.5*(pi - ϕ_1 - ϕ_2)
	θ = pi/2 - γ - ϕ_1

	# length of leg:
	d = WOOFER_CONFIG.l0*sin(γ)
	h1 = WOOFER_CONFIG.l0*cos(γ)
	h2 = sqrt(WOOFER_CONFIG.l1^2 - d^2)
	L = h1 + h2

	unrotated = [L*sin(θ), WOOFER_CONFIG.ABDUCTION_OFFSET, -L*cos(θ)]

	# rotated by abduction
	r_body .= RotX(beta) * unrotated

	if i==1
		r_body .= r_body + [WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
		println(r_body)
	elseif i==2
		r_body .= r_body + [WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
	elseif i==3
		r_body .= r_body + [-WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
	else
		r_body .= r_body + [-WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
	end
end

function legJacobian(α::Vector{Float64})
	ForwardDiff.jacobian(forwardKinematics, α)
end


function forwardKinematicsAll!(r_body::Vector{Float64}, α::Vector{Float64}, abduction_offset=0)
	for i in 1:4
		# leg pointing straight down
		r_i = zeros(3)
		r_i = forwardKinematics!(r_i, α[3*(i-1)+1:3*(i-1)+3], i)

		if i==1
			r_body[1:3] .= r_i + [WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
		elseif i==2
			r_body[4:6] .= r_i + [WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
		elseif i==3
			r_body[7:9] .= r_i + [-WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
		else
			r_body[10:12] .= r_i + [-WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
		end
	end
end

function force2Torque!(τ::Vector{Float64}, f::Vector{Float64}, α::Vector{Float64})
	τ[1:3] .= transpose(legJacobian(α[1:3]))*f[1:3]
	τ[4:6] .= transpose(legJacobian(α[4:6]))*f[4:6]
	τ[7:9] .= transpose(legJacobian(α[7:9]))*f[7:9]
	τ[10:12] .= transpose(legJacobian(α[10:12]))*f[10:12]
end
