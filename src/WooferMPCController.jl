@with_kw mutable struct ControllerParams
	# initialize everything once
   mpc_torques = zeros(12)
   swing_torques = zeros(12)
   swing_torque_i = zeros(3)

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

   N

   mpc_update = 0.001

   contacts = zeros(Int64, 4, N+1)
   foot_locs = zeros(12, N+1)

   x_ref = zeros(12, N+1)

   forces = [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4

   x_des = [0, 0, 0.32, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0, 0.0, 0.00, 0]
end

function mpcControlWoofer!(torques::Vector{T}, x_est::Vector{T}, t::T, joint_pos::Vector{T}, joint_vel::Vector{T}, controller_params::ControllerParams, gait::GaitParams, swing_params::SwingLegParams, footstep_config::FootstepPlannerParams, mpc_config::MPCControllerParams) where {T<:Number}
	# get current leg positions
	forwardKinematicsAll!(controller_params.cur_foot_loc, joint_pos)

	# prev phase -> cur_phase check contacts to regenerate swing
	controller_params.cur_phase = getPhase(t, gait)
	controller_params.active_feet = gait.contact_phases[:, controller_params.cur_phase]
	coordinateExpander!(controller_params.active_feet_12, controller_params.active_feet)

	# swing leg
	for i in 1:4
	   # calculate footstep and generate trajectory (stored in swing params) if needed
	   if gait.contact_phases[i, controller_params.prev_phase] == 1
			if gait.contact_phases[i, controller_params.cur_phase] == 0
	         nextFootstepLocation!(view(swing_params.next_foot_loc, (3*(i-1)+1):(3*(i-1)+3)), controller_params.cur_foot_loc[(3*(i-1)+1):(3*(i-1)+3)], x_est[7:9], x_est[12], gait, nextPhase(controller_params.cur_phase, gait), i)

	         # make sure MPC accounts for this next foot location
	         footstep_config.next_foot_locs[(3*(i-1)+1):(3*(i-1)+3)] .= swing_params.next_foot_loc[(3*(i-1)+1):(3*(i-1)+3)]

	         print("Previous: ")
	         print(controller_params.cur_foot_loc[(3*(i-1)+1):(3*(i-1)+3)])
	         println()
	         print("New: ")
	         print(swing_params.next_foot_loc[(3*(i-1)+1):(3*(i-1)+3)])
	         println()
	         generateFootTrajectory(controller_params.cur_foot_loc[(3*(i-1)+1):(3*(i-1)+3)], x_est[7:9], t, t+gait.phase_times[controller_params.cur_phase], i, swing_params)
	      end
	   end

	   # actually calculate swing torques
	   if gait.contact_phases[i, controller_params.cur_phase] == 0
	      # calculate current foot tip velocity
		  J = zeros(3,3)
		  legJacobian!(J, joint_pos[(3*(i-1)+1):(3*(i-1)+3)])
	      controller_params.cur_foot_vel_i = J * joint_vel[(3*(i-1)+1):(3*(i-1)+3)]

	      calcSwingTorques!(controller_params.swing_torque_i, controller_params.cur_foot_loc[(3*(i-1)+1):(3*(i-1)+3)], controller_params.cur_foot_vel_i, joint_pos[(3*(i-1)+1):(3*(i-1)+3)], t, i, swing_params)
	      controller_params.swing_torques[(3*(i-1)+1):(3*(i-1)+3)] .= controller_params.swing_torque_i
	   end
	end
	controller_params.prev_phase = controller_params.cur_phase

	if t % controller_params.mpc_update < 1e-3
	   # update MPC forces
	   generateReferenceTrajectory!(controller_params.x_ref, x_est, controller_params.x_des, mpc_config)
	   constructFootHistory!(controller_params.contacts, controller_params.foot_locs, t, controller_params.x_ref, controller_params.cur_foot_loc, mpc_config, gait, footstep_config)
	   solveFootForces!(controller_params.forces, controller_params.x_ref, controller_params.contacts, controller_params.foot_locs, mpc_config, true)
	   println(controller_params.forces)
	end

	# needs to be negative so force is exerted by body on world
	force2Torque!(controller_params.mpc_torques, -controller_params.forces, joint_pos)

	torques .= controller_params.active_feet_12 .* controller_params.mpc_torques + (ones(12)-controller_params.active_feet_12) .* controller_params.swing_torques
end
