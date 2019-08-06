using StaticArrays
using Parameters
using LinearAlgebra
using SparseArrays

include("QPBalanceController.jl")
include("WooferDynamics.jl")
include("WooferConfig.jl")
include("ConventionalRodriguesParameters.jl")
include("StandingPlanner.jl")

config = initBalanceController()

t = 0.0
p_ref = zeros(3)
o_ref = zeros(3)

standingPlanner!(p_ref, o_ref, t)

torques = zeros(12)
x = [0, 0, 0.32, 0, 0, 0, 0, 0, 0, 0, 0, 0]

joint_pos = zeros(12)
forces = -[0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4

balanceController!(torques, x, joint_pos, p_ref, o_ref, forces, config)
