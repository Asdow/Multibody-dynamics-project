# File for testing and interactive development
#######################################
#workspace()
include("Datatypes.jl")
include("Body_Init_functions.jl")
include("Solids_functions.jl")
include("Volume_functions.jl")
include("Moment_of_inertia_functions.jl")
include("Coordinate_transform_functions.jl")
include("Quaternion_functions.jl")
include("Plotting_functions.jl")
include("Dynamics_functions.jl")
include("Contact_modeling.jl")

function test()
      # gravity
      g = [0.0; -9.81; 0.0]
      # Init Ground object
      Ground = func_init_body_empty(4, 1, 5.0, Float64)
      fill_rigidbody_Ground!(Ground, 5.0, 5.0, -1.0)
      # Init Array of bodies
      Nbodies = Array(RigidBody, 4)
      Nbodies[1] = Ground
      for i = collect(2:length(Nbodies))
            Nbodies[i] = func_init_body_empty(4, 2, 0.5, Float64)
            fill_rigidbody_cube!(Nbodies[i], 0.5, Nbodies[i].pd.w, 0.5, 2.0)
            # Randomize start location
            Nbodies[i].sv.x[1] = rand(1:50)/100
            Nbodies[i].sv.x[2] = rand(1:200)/100
            # Randomize orientation
            for j in range(1,4)
                  Nbodies[i].sv.q[j] = rand(0:100)/100
            end
            func_euler_normalize!(Nbodies[i].sv.q)
            func_rotation_matrix!(Nbodies[i].aux.Rot_mat,Nbodies[i].sv.q)
            Nbodies[i].aux.Rot_mat_inv = inv(Nbodies[i].aux.Rot_mat)
            # Update world coordinates
            func_body_to_world(Nbodies[i])
            # Update inertias
            func_inertia_body_to_world!(Nbodies[i])
      end
      # init contact list
      contact_list = Array(Contact, 0)
      # init plot
      fig, ax = func_init_plot(1.0, 1.0, 1.0)
      # Draw axes
      origin_axes = func_origin_axes(ax)
      # init array for plotting bodies and draw bodies
      body_plot = Array(PyCall.PyObject, length(Nbodies))
      func_draw_bodies!(ax, body_plot, Nbodies)
      body_axes_array = Array(Array{PyCall.PyObject, 1}, length(Nbodies)-1)
      for i in collect(1:length(Nbodies)-1)
            body_axes_array[i] = func_origin_axes(ax, Nbodies[i+1])
      end

      # Simulation times
      t_init = 0.0
      delta_t = 0.008
      t_end = 8.50
      sleep(0.2)

      # SIMULATION LOOP
      for t = collect(t_init:delta_t:t_end)
            # Clear forces
            Threads.@threads for i in range(1, length(Nbodies))
                  func_clear_forces!(Nbodies[i])
            end
            # Collision handling
            # Check for collisions
            func_ground_contact!(Nbodies, contact_list)
            # Ground collisions
            if length(contact_list) > 0
                  func_ground_collisions!(contact_list, Nbodies)
            end

            # Translation & Rotation
            Threads.@threads for i = collect(1:length(Nbodies))
                  # Translation
                  F_external = g*Nbodies[i].md.massa
                  Nbodies[i].f.force_world += F_external
                  Nbodies[i].sv.P += Nbodies[i].f.force_world*delta_t
                  Nbodies[i].dv.ẋ = Nbodies[i].sv.P*Nbodies[i].md.massa_inv
                  Nbodies[i].sv.x += Nbodies[i].dv.ẋ.*delta_t
                  # Rotation
                  Nbodies[i].sv.L += Nbodies[i].f.torque_world*delta_t
                  Nbodies[i].aux.ω_world = Nbodies[i].J.world_inv*Nbodies[i].sv.L
                  # Aux variables
                  func_omega_tilde4_world!(Nbodies[i].aux.ω_tilde4, Nbodies[i].aux.ω_world)
                  func_euler_vel!(Nbodies[i].dv.q̇, Nbodies[i].sv.q, Nbodies[i].aux.ω_tilde4)
                  func_euler!(Nbodies[i].dv.q̈, Nbodies[i].dv.q̇, Nbodies[i].sv.q, delta_t)
                  func_euler_normalize!(Nbodies[i].sv.q)
                  func_rotation_matrix!(Nbodies[i].aux.Rot_mat,Nbodies[i].sv.q)
                  Nbodies[i].aux.Rot_mat_inv = inv(Nbodies[i].aux.Rot_mat)
                  func_inertia_body_to_world!(Nbodies[i])
                  # Update world coordinates
                  func_body_to_world(Nbodies[i])
            end
            # Remove old plots and update with new ones
            for i = collect(1:length(Nbodies))
                  body_plot[i][:remove]()
                  if i != 1
                        body_axes_array[i-1][1][:remove]()
                        body_axes_array[i-1][2][:remove]()
                        body_axes_array[i-1][3][:remove]()
                  end
            end

            func_draw_bodies!(ax, body_plot, Nbodies)
            for i in collect(1:length(Nbodies)-1)
                  func_origin_axes!(ax, Nbodies[i+1], body_axes_array[i])
            end

            fig[:canvas][:flush_events]()
      end
      return nothing
end
