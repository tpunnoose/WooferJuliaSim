using ForwardDiff

# function forwardKinematics(α::Vector{T}) where {T<:Real}
# 	beta = α[1]
# 	ϕ_1 = α[2]
# 	ϕ_2 = α[3]
#
#
# 	# coordinates of leg pointing straight down
# 	γ = 0.5*(pi - ϕ_1 + ϕ_2)
# 	θ = pi/2 - γ - ϕ_1
#
# 	# length of leg:
# 	d = WOOFER_CONFIG.l0*sin(γ)
# 	h1 = WOOFER_CONFIG.l0*cos(γ)
# 	h2 = sqrt(WOOFER_CONFIG.l1^2 - d^2)
# 	L = h1 + h2
#
# 	unrotated = [L*sin(θ), WOOFER_CONFIG.ABDUCTION_OFFSET, -L*cos(θ)]
#
# 	# rotated by abduction
# 	return RotX(beta) * unrotated
# end

function forwardKinematics!(r_body::Vector{T}, α::Vector{T}, i::Int64=1) where {T<:Real}
	beta = α[1]
	ϕ_1 = α[2] # TODO: check the sign of these
	ϕ_2 = α[3] # TODO: check the sign of these


	# coordinates of leg pointing straight down
	γ = 0.5*(pi - ϕ_1 + ϕ_2)
	θ = pi/2 - γ - ϕ_1

	# length of leg:
	d = WOOFER_CONFIG.l0*sin(γ)
	h1 = WOOFER_CONFIG.l0*cos(γ)
	h2 = sqrt(WOOFER_CONFIG.l1^2 - d^2)
	L = h1 + h2

	if i==1
		unrotated = [L*sin(θ), -WOOFER_CONFIG.ABDUCTION_OFFSET, -L*cos(θ)]
		r_body .= RotX(beta) * unrotated + [WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
	elseif i==2
		unrotated = [L*sin(θ), WOOFER_CONFIG.ABDUCTION_OFFSET, -L*cos(θ)]
		r_body .= RotX(beta) * unrotated + [WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
	elseif i==3
		unrotated = [L*sin(θ), -WOOFER_CONFIG.ABDUCTION_OFFSET, -L*cos(θ)]
		r_body .= RotX(beta) * unrotated + [-WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
	else
		unrotated = [L*sin(θ), WOOFER_CONFIG.ABDUCTION_OFFSET, -L*cos(θ)]
		r_body .= RotX(beta) * unrotated + [-WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
	end
end

# function legJacobian(α::Vector{Float64})
# 	ForwardDiff.jacobian(forwardKinematics, α)
# end

function legJacobian!(jacobian::Array{T,2}, α::Vector{T}) where {T<:Real}
	r_body = zeros(3)
	ForwardDiff.jacobian!(jacobian, forwardKinematics!, r_body, α)
end


function forwardKinematicsAll!(r_body::Vector{Float64}, α::Vector{Float64})
	r_i = zeros(3)
	for i in 1:4
		forwardKinematics!(r_i, α[3*(i-1)+1:3*(i-1)+3], i)
		r_body[3*(i-1)+1:3*(i-1)+3] .= r_i
	end
end

function force2Torque!(τ::Vector{Float64}, f::Vector{Float64}, α::Vector{Float64})
	J = zeros(3,3)
	for i in 1:3:12
		legJacobian!(J, α[i:i+2])
		τ[i:i+2] .= transpose(J)*f[i:i+2]
	end
end
