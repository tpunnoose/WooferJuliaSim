using LinearAlgebra
using Parameters
using Rotations

include("WooferDynamics.jl")
include("WooferConfig.jl")

x_dot = zeros(12)
x = [0.0, 0.0, 0.32, zeros(9)...]
u = [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4
joint_pos = zeros(12)

A = zeros(12,12)
B = zeros(12,12)

A!(A, x, u, joint_pos)
# B!(B, x, u, joint_pos)
