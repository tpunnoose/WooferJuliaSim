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
   x_true = zeros(12)
   x_est = zeros(12)
   x = zeros(13)
   g = zeros(3)
   accel = zeros(3)
   gyro = zeros(3)
   joint_pos = zeros(12)
   joint_vel = zeros(12)
   contacts = zeros(4)

   ψ = 0.0
   # x_des = [0.05, 0.0, 0.32, 0.05, -0.02, ψ, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0]
   x_des = [0.0, 0.0, 0.32, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0]

   lower_dt = 0.001
   # x0 = [0.34, zeros(11)...]
   # P = Diagonal([0.01, 0.1*ones(2)..., 0.1*ones(3)..., 0.1*ones(6)...])
   # Q = Diagonal([0.001, 0.001*ones(2)..., 0.001*ones(3)..., 1e-10*ones(6)...])
   # R = Diagonal([0.001, 0.01*ones(3)...])
   # est_params = StateEstimatorParams(dt=lower_dt, x=x0, P=P, Q=Q, R=R)
   # last_t = 0.0

   qp_params = initQPParams(lower_dt, x_des, ψ)
   qp_forces = zeros(12)
   qp_torques = zeros(12)

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

               q = Quat(x[4], x[5], x[6], x[7])
               R = SMatrix{3,3}(q)

               x_b = R'*x[1:3]
               v_b = R'*x[8:10]

               ## Add in state estimation here ##
               x_true[1:3] .= s.d.qpos[1:3]
               x_true[4:6] .= 2*s.d.qpos[5:7] # is this right q -> g?
               x_true[7:9] .= s.d.qvel[1:3]
               x_true[10:12] .= s.d.qvel[4:6] # is this in the body frame?

               accel       .= s.d.sensordata[1:3]
               gyro        .= s.d.sensordata[4:6]
               joint_pos   .= s.d.sensordata[7:18]
               joint_vel   .= s.d.sensordata[19:30]

               # stateEstimatorUpdate(t-last_t, accel, gyro, joint_pos, joint_vel, contacts, est_params)
               # x_est[1:6] .= est_params.x[1:6]
               # x_est[7:9] .= gyro
               # last_t = t

               # QP Balance Controller
               qpBalance!(qp_forces, x_true, qp_params)
               @time force2Torque!(qp_torques, -(qp_forces), joint_pos)
               s.d.ctrl .= qp_torques
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
