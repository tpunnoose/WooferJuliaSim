function simulate()
   # look at local directory
   pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__)

   # parse MuJoCo XML file
   xmlparser = pyimport("WooferXMLParser")
   xmlparser.Parse()

   s = loadmodel("src/woofer_out.xml", 1200, 900)

   d = s.d
   m = s.m

   # initialize everything once
   mpc_torques = zeros(12)
   x_true = zeros(9)
   x_est = zeros(9)
   x = zeros(13)
   p_ref = zeros(3)
   o_ref = zeros(3)
   g = zeros(3)
   accel = zeros(3)
   gyro = zeros(3)
   joint_pos = zeros(12)
   joint_vel = zeros(12)
   contacts = zeros(4)

   # planning_dt = 0.04
   # N = 15
   # mpc_update = 0.01
   planning_dt = 0.001
   N = 1
   mpc_update = 0.001

   x_des = [0.32, 0.00, 0.0, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0]

   lower_dt = 0.001
   x0 = [0.34, zeros(11)...]
   P = Diagonal([0.01, 0.1*ones(2)..., 0.1*ones(3)..., 0.1*ones(6)...])
   Q = Diagonal([0.001, 0.001*ones(2)..., 0.001*ones(3)..., 1e-10*ones(6)...])
   R = Diagonal([0.001, 0.01*ones(3)...])
   est_params = StateEstimatorParams(dt=lower_dt, x=x0, P=P, Q=Q, R=R)
   last_t = 0.0

   qp_config = initBalanceController()
   qp_torques = zeros(12)
   qp_forces = -[0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0]*WOOFER_CONFIG.MASS*9.81/4
   x_qp = zeros(12)

   lqr_params = initLQRParams(lower_dt, x_des[1:9])
   lqr_forces = zeros(12)
   lqr_torques = zeros(12)

   plot_length = 10000
   plot_est_state = zeros(12, plot_length)
   plot_true_state = zeros(6, plot_length)
   plot_i = 1

   mpc_config = initMPCControllerConfig(planning_dt, N)
   # gait = GaitParams(num_phases=1, contact_phases=[1;1;1;1], phase_times=[1.0]) # standing gait
   gait = GaitParams(num_phases=1, contact_phases=[1;0;1;0], phase_times=[1.0]) # standing 3 legged gait
   # gait = GaitParams() # trot gait
   footstep_config = FootstepPlannerParams()
   swing_params = SwingLegParams()
   controller_params = ControllerParams(N=N, mpc_update=mpc_update, x_des=x_des)

   # Loop until the user closes the window
   WooferMocap.alignscale(s)
   while !GLFW.WindowShouldClose(s.window)
      ### basically sim step so things don't have to be defined multiple times
      if s.paused
         if s.pert[].active > 0
            mjv_applyPerturbPose(m, d, s.pert, 1)  # move mocap and dynamic bodies
            mj_forward(m, d)
         end
      else
         #slow motion factor: 10x
         factor = s.slowmotion ? 10 : 1

         # advance effective simulation time by 1/refreshrate
         startsimtm = d.d[].time
         starttm = time()
         refreshtm = 1.0/(factor*s.refreshrate)
         updates = refreshtm / m.m[].opt.timestep

         steps = round(Int, round(s.framecount+updates)-s.framecount)
         s.framecount += updates

         for i=1:steps
            # clear old perturbations, apply new
            d.xfrc_applied .= 0.0

            # add in noise like perturbations
            # d.xfrc_applied[7:9] .= [5, 5, -10]
            # d.xfrc_applied[7:9] .= 10*randn(Float64, 3)
            # d.xfrc_applied[10:12] .= 5*randn(Float64, 3)

            if s.pert[].select > 0
               mjv_applyPerturbPose(m, d, s.pert, 0) # move mocap bodies only
               mjv_applyPerturbForce(m, d, s.pert)
            end

            t = d.d[].time

            # lower level update loop (eg state estimation, torque updates)
            if t % lower_dt < 1e-3
               # ground truth states
               x[1:3] .= s.d.qpos[1:3]
               x[4:7] .= s.d.qpos[4:7]
               x[8:10] .= s.d.qvel[1:3]
               x[11:13] .= s.d.qvel[4:6]

               q = Quat(x[4], x[5], x[6], x[7])
               R = SMatrix{3,3}(q)

               x_b = R'*x[1:3]
               v_b = R'*x[8:10]

               ## Add in state estimation here ##
               x_true[1] = x_b[3]
               x_true[2] = atan(2*(s.d.qpos[4]*s.d.qpos[5] + s.d.qpos[6]*s.d.qpos[7])/(1 - 2*(s.d.qpos[5]^2 + s.d.qpos[6]^2)))
               x_true[3] = asin(2*(s.d.qpos[4]*s.d.qpos[6] - s.d.qpos[5]*s.d.qpos[7]))
               x_true[4:6] .= v_b
               x_true[7:9] .= x[11:13]

               accel       .= s.d.sensordata[1:3]
               gyro        .= s.d.sensordata[4:6]
               joint_pos   .= s.d.sensordata[7:18]
               joint_vel   .= s.d.sensordata[19:30]
               # wait for feet to actually be in contact
               if t > 0.1
                  contacts .= gait.contact_phases[:,getPhase(t, gait)]
               else
                  contacts .= zeros(4)
               end

               stateEstimatorUpdate(t-last_t, accel, gyro, joint_pos, joint_vel, contacts, est_params)
               x_est[1:6] .= est_params.x[1:6]
               x_est[7:9] .= gyro
               last_t = t

               if plot_i < plot_length
                  # plot_est_state[1:9, plot_i] .= x_est[1:9]
                  plot_est_state[:, plot_i] .= est_params.x
                  plot_true_state[:, plot_i] .= x_true[1:6]
                  plot_i += 1
               end

               # # MPC Controller
               # mpcControlWoofer!(mpc_torques, x_true, t, joint_pos, joint_vel, controller_params, gait, swing_params, footstep_config, mpc_config)
               # s.d.ctrl .= mpc_torques

               # QP Controller
               # standingPlanner!(p_ref, o_ref, t)
               # quaternion2CRP!(g, x[4:7])
               # x_qp[1:3] .= x[1:3]
               # x_qp[4:6] .= g
               # x_qp[7:9] .= x[8:10]
               # x_qp[10:12] .= x[11:13]
               # balanceController!(qp_torques, x_qp, joint_pos, p_ref, o_ref, qp_forces, qp_config)
               # s.d.ctrl .= qp_torques

               # LQR Balance Controller
               # lqrBalance!(lqr_forces, x_est, joint_pos, lqr_params)
               lqrBalance!(lqr_forces, x_true, joint_pos, lqr_params)
               force2Torque!(lqr_torques, -lqr_forces, joint_pos)
               s.d.ctrl .= lqr_torques
            end

            mj_step(s.m, s.d)

            # break on reset
            (d.d[].time < startsimtm) && break
         end
      end

      render(s, s.window)
      GLFW.PollEvents()
   end
   GLFW.DestroyWindow(s.window)

   print("RMS roll error: ")
   println(sqrt(mean((plot_est_state[2,:] - plot_true_state[2,:]).^2)))

   print("RMS pitch error: ")
   println(sqrt(mean((plot_est_state[3,:] - plot_true_state[3,:]).^2)))

   print("RMS vx error: ")
   println(sqrt(mean((plot_est_state[4,:] - plot_true_state[4,:]).^2)))

   print("RMS vy error: ")
   println(sqrt(mean((plot_est_state[5,:] - plot_true_state[5,:]).^2)))

   print("RMS vz error: ")
   println(sqrt(mean((plot_est_state[6,:] - plot_true_state[6,:]).^2)))

   # plot(plot_est_state[2,:])
   # plot!(plot_true_state[2,:])
   #
   # plot(plot_est_state[3,:])
   # plot!(plot_true_state[3,:])

   # plot(plot_est_state[4,:])
   # plot!(plot_true_state[4,:])

   # plot(plot_est_state[5,:])
   # plot!(plot_true_state[5,:])
   #
   # plot(plot_est_state[6,:])
   # plot!(plot_true_state[6,:])
end
