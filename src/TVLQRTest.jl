using LinearAlgebra
using Parameters
using Rotations

include("WooferDynamics.jl")
include("WooferConfig.jl")

x_dot = zeros(12)
x = [0.0, 0.0, 0.32, zeros(9)...]
u = zeros(12)#[0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4
joint_pos = zeros(12)

f = nonlinearDynamics([x..., u..., joint_pos...])

ForwardDiff.jacobian(nonlinearDynamics, [x..., u..., joint_pos...])

(A, B) = AB(x, u, joint_pos)

A

B
