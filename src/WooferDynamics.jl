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

"""
create Rotation matrix from quaternion projection
"""
function rotationMatrix!(Q::Array{T,2}, x::Vector{T}) where {T<:Real}

end

"""
Full nonlinear dynamics for the rigid body model of Woofer
"""
function nonlinearDynamics!(x_dot::Vector{T}, x::Vector{T}, u::Vector{T}, joint_pos::Vector{T}) where {T<:Real}
	x_dot .= zeros(12)
	x_dot[1:3] .= x[7:9]

	# orientation dynamics
	x_dot[4] = x[10]*sqrt(1-0.25*(x[4]^2 + x[5]^2 + x[6]^2)) + 0.5*(x[5]*x[12] - x[6]*x[11])
	x_dot[5] = x[11]*sqrt(1-0.25*(x[4]^2 + x[5]^2 + x[6]^2)) + 0.5*(x[4]*x[12] - x[6]*x[10])
	x_dot[6] = x[12]*sqrt(1-0.25*(x[4]^2 + x[5]^2 + x[6]^2)) + 0.5*(x[4]*x[11] - x[5]*x[10])

	# TODO: take out allocations here
	r_i = zeros(3)
	r_hat = zeros(3,3)
	ω_hat = zeros(3,3)
	Q = Matrix{Float64}(I, 3, 3)

	# TODO: add rotation matrix here as function of x[4:6]
	# rotationMatrix!(Q, x[4:6])

	x_dot[7:9] .= -[0, 0, 9.81]
	x_dot[10:12] .= -WOOFER_CONFIG.INV_INERTIA*ω_hat*WOOFER_CONFIG.INERTIA*x[10:12]

	for i=1:4
		# forces map to accelerations
		x_dot[7:9] .= x_dot[7:9] + 1/WOOFER_CONFIG.MASS*u[3*(i-1)+1:3*(i-1)+3]

		# torques map to angular accelerations
		forwardKinematics!(r_i, joint_pos[i:i+2], i)
		skewSymmetricMatrix!(r_hat, r_i)
		x_dot[10:12] .= x_dot[10:12] + WOOFER_CONFIG.INV_INERTIA*r_hat*Q*u[3*(i-1)+1:3*(i-1)+3]
	end
end

"""
Not inplace version of nonlinearDynamics!

z = [x u joint_pos]
"""
function nonlinearDynamics(z::Vector{T}) where {T<:Real}
	x_dot = zeros(12)

	###

	x_dot[1:3] .= x[7:9]

	# orientation dynamics
	x_dot[4] = x[10]*sqrt(1-0.25*(x[4]^2 + x[5]^2 + x[6]^2)) + 0.5*(x[5]*x[12] - x[6]*x[11])
	x_dot[5] = x[11]*sqrt(1-0.25*(x[4]^2 + x[5]^2 + x[6]^2)) + 0.5*(x[4]*x[12] - x[6]*x[10])
	x_dot[6] = x[12]*sqrt(1-0.25*(x[4]^2 + x[5]^2 + x[6]^2)) + 0.5*(x[4]*x[11] - x[5]*x[10])

	# TODO: take out allocations here
	r_i = zeros(3)
	r_hat = zeros(3,3)
	ω_hat = zeros(3,3)
	Q = Matrix{Float64}(I, 3, 3)

	# TODO: add rotation matrix here as function of x[4:6]
	# rotationMatrix!(Q, x[4:6])

	x_dot[7:9] .= -[0, 0, 9.81]
	x_dot[10:12] .= -WOOFER_CONFIG.INV_INERTIA*ω_hat*WOOFER_CONFIG.INERTIA*x[10:12]

	for i=1:4
		# forces map to accelerations
		x_dot[7:9] .= x_dot[7:9] + 1/WOOFER_CONFIG.MASS*u[3*(i-1)+1:3*(i-1)+3]

		# torques map to angular accelerations
		forwardKinematics!(r_i, joint_pos[i:i+2], i)
		skewSymmetricMatrix!(r_hat, r_i)
		x_dot[10:12] .= x_dot[10:12] + WOOFER_CONFIG.INV_INERTIA*r_hat*Q*u[3*(i-1)+1:3*(i-1)+3]
	end

	###
	return x_dot
end

"""
Create linearized A and B matrix
"""
function AB(x::Vector, u::Vector, joint_pos::Vector)
	J = ForwardDiff.jacobian(nonlinearDynamics, [x..., u..., joint_pos...])

	A = J[1:12,1:12]
	B = J[1:12,13:24]
	return (A, B)
end

"""
Fills A with skew symmetric version of a
"""
function skewSymmetricMatrix!(A::Matrix, a::Vector)
	A[1,2] = -a[3]
	A[1,3] = a[2]
	A[2,1] = a[3]
	A[2,3] = -a[1]
	A[3,1] = -a[2]
	A[3,2] = a[1]
end
