using Revise

include("./src/WooferSim.jl")

using .WooferSim
WooferSim.simulate()
