# Woofer Julia Simulation

## Overview
This repository contains the most up to date Woofer code that accounts for full state feedback from the MOCAP room.

## Installation for Simulation
1. Acquire a license for MuJoCo at http://mujoco.org/. You can get a free trial of the professional license for a month, or with a student account, a free year.

2. Save the license ```mjkey.txt``` somewhere and set the environment variable ```MUJOCO_KEY_PATH``` to that location. One way to set the environment variable is through your bash profile. On a mac this is done by adding the line
```
export MUJOCO_KEY_PATH=[YOUR PATH]/mjkey.txt
```
to your ~/.bash_profile.

3. Install Julia by visiting https://julialang.org/downloads/

4. Install MuJoCo.jl. Enter the Package Manager by typing the `]` key and type  
```shell
pkg> add https://github.com/klowrey/MuJoCo.jl
```
5. Install this repository through Pkg
```shell
pkg> add https://github.com/tpunnoose/WooferSim
```

## Run Simulation
1. Enter the Julia REPL in the WooferSim directory (this will be in the .julia/dev folder).
2. Run
```julia
include("main.jl")
```
In order to change what simulation is being run, scroll to the bottom of the WooferSim.jl and uncomment the desired controller.
3. The MuJoCo simulator should then pop up in a new window with various interactive options. Press space to start the simulation.
- Click and drag with the left mouse button to orbit the camera, and with the right mouse button to pan the camera.
- To perturb the robot, double click on the body you want to perturb, then hold Control and click and drag with the mouse. Using the left mouse button will apply a rotational torque while the right button will apply a translational force.
- Press space to play or pause the simulation.
