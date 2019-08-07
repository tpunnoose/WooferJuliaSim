using StaticArrays
using Parameters
using LinearAlgebra
using SparseArrays
using Rotations

include("WooferConfig.jl")
include("MPCController.jl")
include("WooferDynamics.jl")
include("Gait.jl")
include("FootstepPlanner.jl")

dt = 0.01
N = 15

mpc_config = initMPCControllerConfig(dt, N, WOOFER_CONFIG)
standingGait = GaitParams(num_phases=1, contact_phases=[1;1;1;1], phase_times=[1.0])
footstep_config = FootstepPlannerParams()

x0 = [0.3, 0.01, 0.01, 0.01, -0.05, 0.00, 0.001, -0.001, 0.04, 9.81]
x_des = [0.31, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0, 0.0, 0.00, 9.81]

x_ref = zeros(10, N)

generateReferenceTrajectory!(x_ref, x0, x_des, mpc_config)

contacts = zeros(Int64, 4, N)
foot_locs = zeros(12, N)

r1 = zeros(3)
r2 = zeros(3)
r3 = zeros(3)
r4 = zeros(3)

forwardKinematics!(r1, zeros(3), 1)
forwardKinematics!(r2, zeros(3), 2)
forwardKinematics!(r3, zeros(3), 3)
forwardKinematics!(r4, zeros(3), 4)

cur_foot_loc = zeros(12)
cur_foot_loc[1:3] = r1
cur_foot_loc[4:6] = r2
cur_foot_loc[7:9] = r3
cur_foot_loc[10:12] = r4

t = 0.0

constructFootHistory!(contacts, foot_locs, t, x_ref, cur_foot_loc, mpc_config, standingGait, footstep_config)

foot_locs

forces = zeros(12*mpc_config.N)

@time solveFootForces!(forces, x0, x_ref, contacts, foot_locs, mpc_config, WOOFER_CONFIG)

forces[1:12]

x_exp = mpc_config.A_qp*x0 + mpc_config.B_qp*forces

x_exp[141:150]
