using StaticArrays
using Parameters
using LinearAlgebra
using SparseArrays

include("WooferConfig.jl")
include("MPCController.jl")
include("WooferDynamics.jl")

dt = 0.01
N = 15

mpc_config = initMPCControllerConfig(dt, N, WOOFER_CONFIG)

x0 = [0.3, 0.01, 0.01, 0.02, -0.01, 0.01, 0.001, -0.001, 0.02, 9.81]

x_ref = repeat(x0, 1, mpc_config.N)
contacts = repeat([1, 1, 1, 1], 1, mpc_config.N)

r_1 =
foot_locs = repeat()
