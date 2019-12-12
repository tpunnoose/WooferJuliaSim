using Parameters
using LinearAlgebra

include("../src/WooferDynamics.jl")
include("../src/WooferConfig.jl")

@with_kw mutable struct StateEstimatorParams
    x::Vector{Float64}
    P::Matrix{Float64}

    xdot::Vector{Float64} = zeros(size(x)[1])
    x_::Vector{Float64} = zeros(size(x)[1])
    x_plus::Vector{Float64} = zeros(size(x)[1])
    P_::Matrix{Float64} = zeros(size(x)[1], size(x)[1])
    P_plus::Matrix{Float64} = zeros(size(x)[1], size(x)[1])

    A::Matrix{Float64} = zeros(size(x)[1], size(x)[1])
    C::Matrix{Float64} = zeros(4, size(x)[1])
    S::Matrix{Float64} = zeros(4, 4)
    y_meas::Vector{Float64} = zeros(4)
    y_pred::Vector{Float64} = zeros(4)
    nu::Vector{Float64} = zeros(4)
    K::Matrix{Float64} = zeros(size(x)[1], 4)
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

    r_fr::Vector{Float64} = [WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
    r_fl::Vector{Float64} = [WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
    r_br::Vector{Float64} = [-WOOFER_CONFIG.LEG_FB, -WOOFER_CONFIG.LEG_LR, 0]
    r_bl::Vector{Float64} = [-WOOFER_CONFIG.LEG_FB, WOOFER_CONFIG.LEG_LR, 0]
end

function SkewSymmetricMatrix(q::Vector{T}) where {T<:Real}
    return [0.0 -q[3] q[2]; q[3] 0.0 -q[1]; -q[2] q[1] 0.0]
end

function L_q(q::Vector{T}) where {T<:Real}
    return [q[1] -q[2:4]'; q[2:4] q[1] * I + SkewSymmetricMatrix(q[2:4])]
end

function R_q(q::Vector{T}) where {T<:Real}
    return [q[1] -q[2:4]'; q[2:4] q[1] * I - SkewSymmetricMatrix(q[2:4])]
end

function InitStateEstimator(
    dt::T,
    x::Vector{T},
    P::Diagonal{T},
    Q::Diagonal{T},
    R::Diagonal{T},
) where {T<:Number} # default Q, R?
    return StateEstimatorParams(dt = dt, x = x, P = P, Q = Q, R = R)
end

function StateEstimatorUpdate(
    dt::AbstractFloat,
    r_hat::Vector,
    q::Vector,
    joint_pos::Vector,
    joint_vel::Vector,
    est_params::StateEstimatorParams,
    lqr_forces::Vector,
)
    x_r = est_params.x[1:3]
    x_ϕ = est_params.x[4:6]
    x_v = est_params.x[7:9]
    x_w = est_params.x[10:12]

    # Reconstruct the orientation quaternion from our state estimate
    sc_mag = 1.0 - x_ϕ[1]^2 - x_ϕ[2]^2 - x_ϕ[3]^2
    x_ϕ_sc = sign(sc_mag) * sqrt(abs(sc_mag))

    # Orientation from mocap
    q_sc = q[1]
    q_hat = q[2:4]

    # Mass and inertia of the robot body
    m = WOOFER_CONFIG.MASS
    J = WOOFER_CONFIG.INERTIA

    # Gives you last 3 components of the quaternion
    V = zeros(3, 4)
    V[:, 2:4] = Matrix{Float64}(I, 3, 3)

    # A is continuous, xdot = Ax
    # Note that we will take gravity into account in the calculation of the net force on the body
    A = zeros(12, 12)
    A[1:3, 7:9] = Matrix{Float64}(I, 3, 3)

    # How qdot depends on q (and omega)
    A[4:6, 4:6] = 0.5 * V * R_q([0; x_w]) * [x_ϕ' / sqrt(abs(1 - x_ϕ' * x_ϕ)); I]

    # How qdot depends on qdot
    A[4:6, 10:12] = 0.5 * (L_q([x_ϕ_sc; x_ϕ])*V')[2:4, :]

    # How qdd depends on qdot
    A[10:12, 10:12] = inv(J) * (SkewSymmetricMatrix(x_w) * J)

    # Control mapping matrix
    B = zeros(12, 6)
    B[7:9, 1:3] = Matrix{Float64}(I, 3, 3) / m
    B[10:12, 4:6] = inv(J)

    # Convert from continuous system to discrete system
    cont_sys = zeros(18, 18)
    cont_sys[1:12, 1:12] .= A
    cont_sys[1:12, 13:18] .= B
    disc_sys = exp(cont_sys * dt)
    A_discrete = disc_sys[1:12, 1:12]
    B_discrete = disc_sys[1:12, 13:18]

    # Mapping from state to output.
    C = [I zeros(3, 9); I zeros(3, 9)]

    # Calculate the measurement error
    y_hat = C * est_params.x[1:12]
    y_true = [r_hat; q_hat]

    # Vector from CoM to foot in body frame
    leg_vec = zeros(3)

    jacobian = zeros(3, 3)
    body_F = zeros(3)
    body_τ = zeros(3)
    for i = 1:4
        leg_joint_pos = joint_pos[3*(i-1)+1:3*(i-1)+3]
        legJacobian!(jacobian, leg_joint_pos, i)
        forwardKinematics!(leg_vec, leg_joint_pos, i)

        lqr_forces_leg = lqr_forces[3*(i-1)+1:3*(i-1)+3]
        τ_leg = cross(leg_vec, lqr_forces_leg)
        body_F += lqr_forces_leg
        body_τ += τ_leg
    end

    # Subtract the weight of the robot from the net forces applied by the legs
    body_F -= est_params.g_n * m

    # Combine the force and torque on the body into one vector
    u = [body_F; body_τ]

    # Propogate dynamics and covariance
    # x_ is the a priori estimate of the true state, x is the a posteriori estimate of the true state
    # P_ is the a priori estimate of the state covariance, P is the a posteriori estimate of the true state covariance
    est_params.x_ = A_discrete * est_params.x + B_discrete * u
    est_params.P_ = A_discrete * est_params.P * A_discrete' + est_params.Q

    K_s = est_params.P_ * C' * inv(C * est_params.P_ * C' + est_params.R)
    est_params.x = est_params.x_ + K_s * (y_true - y_hat)
    est_params.P = est_params.P_ - K_s * C * est_params.P_
end
