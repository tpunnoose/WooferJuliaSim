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
   x_true = zeros(12)
   x = zeros(13)
   accel = zeros(3)
   gyro = zeros(3)
   joint_pos = zeros(12)
   joint_vel = zeros(12)
   contacts = zeros(4)

   # planning_dt = 0.04
   # N = 15
   planning_dt = 0.001
   N = 1
   mpc_update = 0.001

   x_des = [0.0, 0.0, 0.32, zeros(9)...]
   # x_des = [0.0, 0.0, 0.4, zeros(3)..., 0.2, 0.0, zeros(4)...]

   lower_dt = 0.001

   mpc_config = initMPCControllerConfig(planning_dt, N)
   gait = GaitParams(num_phases=1, contact_phases=[1;1;1;1], phase_times=[1.0]) # standing gait
   # gait = GaitParams(num_phases=1, contact_phases=[1;0;1;0], phase_times=[1.0]) # standing 3 legged gait
   # gait = GaitParams() # trot gait
   footstep_config = FootstepPlannerParams()
   swing_params = SwingLegParams()
   controller_params = ControllerParams(N=N, mpc_update=mpc_update, x_des=x_des)

   # Loop until the user closes the window
   WooferSim.alignscale(s)
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

               x_true[1:3] .= s.d.qpos[1:3]
               x_true[4:6] .= s.d.qpos[5:7] # is this right q -> g?
               x_true[7:9] .= s.d.qvel[1:3]
               x_true[10:12] .= s.d.qvel[4:6] # is this in the body frame?

               accel       .= s.d.sensordata[1:3]
               gyro        .= s.d.sensordata[4:6]
               joint_pos   .= s.d.sensordata[7:18]
               joint_vel   .= s.d.sensordata[19:30]

               # MPC Controller
               mpcControlWoofer!(mpc_torques, x_true, t, joint_pos, joint_vel, controller_params, gait, swing_params, footstep_config, mpc_config)
               s.d.ctrl .= mpc_torques
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
end
