@with_kw struct ControllerParams
	# initialize everything once
   mpc_torques = zeros(12)
   swing_torques = zeros(12)
   swing_torque_i = zeros(3)
   x_est = zeros(10)
   x_est[10] = 9.81
   x = zeros(13)

   prev_phase = 1
   cur_phase = 1

   r1 = zeros(3)
   r2 = zeros(3)
   r3 = zeros(3)
   r4 = zeros(3)
   cur_foot_loc = zeros(12)
   cur_foot_vel_i = zeros(3)
   active_feet = zeros(Int64, 4)
   active_feet_12 = zeros(Int64, 12)

   dt = 0.01
   N = 10

   mpc_update = 0.01

   mpc_config = initMPCControllerConfig(dt, N, WOOFER_CONFIG)
   # gait = GaitParams(num_phases=1, contact_phases=[1;1;1;1], phase_times=[1.0]) # standing gait
   gait = GaitParams() # trot gait
   footstep_config = FootstepPlannerParams()
   swing_params = SwingLegParams()

   contacts = zeros(Int64, 4, N)
   foot_locs = zeros(12, N)

   x_ref = zeros(10, N)

   forces = -[0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4

   x_des = [0.31, 0.00, 0.0, 0.0, 0.00, 0.00, 0.0, 0.0, 0.0, 9.81]
end

function controlWoofer!(torques::Vector{T}, x_est::Vector{T}, t::T, joint_pos::Vector{T}, joint_vel::Vector{T}, controller_params::ControllerParams) where {T<:Number}
	# get current leg positions
	# TODO: put this in a function
	forwardKinematics!(controller_params.r1, joint_pos[1:3], 1)
	forwardKinematics!(controller_params.r2, joint_pos[4:6], 2)
	forwardKinematics!(controller_params.r3, joint_pos[7:9], 3)
	forwardKinematics!(controller_params.r4, joint_pos[10:12], 4)

	controller_params.cur_foot_loc[1:3] = controller_params.r1
	controller_params.cur_foot_loc[4:6] = controller_params.r2
	controller_params.cur_foot_loc[7:9] = controller_params.r3
	controller_params.cur_foot_loc[10:12] = controller_params.r4

	# prev phase -> cur_phase check contacts to regenerate swing
	controller_params.cur_phase = getPhase(t, controller_params.gait)
	controller_params.active_feet = controller_params.gait.contact_phases[:, controller_params.cur_phase]
	coordinateExpander!(controller_params.active_feet_12, controller_params.active_feet)

	# swing leg
	for i in 1:4
	   # calculate footstep and generate trajectory (stored in swing params) if needed
	   if gait.contact_phases[i, prev_phase] == 1
			if gait.contact_phases[i, cur_phase] == 0
			 nextFootstepLocation!(view(swing_params.next_foot_loc, (3*(i-1)+1):(3*(i-1)+3)), cur_foot_loc[(3*(i-1)+1):(3*(i-1)+3)], x_est[4:6], x_est[9], gait, nextPhase(cur_phase, gait), i)

			 # make sure MPC accounts for this next foot location
			 footstep_config.next_foot_locs[(3*(i-1)+1):(3*(i-1)+3)] .= swing_params.next_foot_loc[(3*(i-1)+1):(3*(i-1)+3)]

			 print("Previous: ")
			 print(cur_foot_loc[(3*(i-1)+1):(3*(i-1)+3)])
			 println()
			 print("New: ")
			 print(swing_params.next_foot_loc[(3*(i-1)+1):(3*(i-1)+3)])
			 println()
			 generateFootTrajectory(cur_foot_loc[(3*(i-1)+1):(3*(i-1)+3)], x_est[4:6], t, t+gait.phase_times[cur_phase], i, swing_params)
		  end
	   end

	   # actually calculate swing torques
	   if gait.contact_phases[i, cur_phase] == 0
		  # calculate current foot tip velocity
		  cur_foot_vel_i = legJacobian(joint_pos[(3*(i-1)+1):(3*(i-1)+3)]) * joint_vel[(3*(i-1)+1):(3*(i-1)+3)]

		  calcSwingTorques!(swing_torque_i, cur_foot_loc[(3*(i-1)+1):(3*(i-1)+3)], cur_foot_vel_i, joint_pos[(3*(i-1)+1):(3*(i-1)+3)], t, i, swing_params)
		  swing_torques[(3*(i-1)+1):(3*(i-1)+3)] .= swing_torque_i
	   end
	end
	prev_phase = cur_phase

	if t % mpc_update < 1e-3
	   # update MPC forces
	   generateReferenceTrajectory!(x_ref, x_est, x_des, mpc_config)
	   constructFootHistory!(contacts, foot_locs, t, x_ref, cur_foot_loc, mpc_config, gait, footstep_config)
	   solveFootForces!(forces, x_est, x_ref, contacts, foot_locs, mpc_config, WOOFER_CONFIG)
	   # println(mpc_torques)
	   # println(forces)
	end

	# needs to be negative so force is exerted by body on world
	force2Torque!(mpc_torques, -forces, joint_pos)

	torques .= active_feet_12 .* mpc_torques + (ones(12)-active_feet_12) .* swing_torques
end
