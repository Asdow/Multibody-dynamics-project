# File for testing and interactive development
#######################################
include("Multiprocess.jl")
include("Packages.jl")
include("Globals.jl")
include("apufunktiot.jl")
include("Quaternions.jl")
include("Datatypes.jl")
include("Visualization.jl")
include("Body_Init_functions.jl")
include("Moment_of_inertia_functions.jl")
include("Solids_functions.jl")
include("Collision_detection.jl")
include("Coordinate_transformations.jl")
include("Dynamics_functions.jl")

function test()
      body = init_body_empty(4, 2, 15.0, 0.6, 0.4);
      # Define body as a cube
      create_cube!(body, 2.0, 2.0, 2.0)
      # set location and orientation to random
      body.sv.x[:] = rand(3)
      body.sv.q[:] = rand(4)
      normalize!(body.sv.q)
      rotmat!(body)
      # Set world coords
      body2world!(body)
      # Create plane so body doesn't fall into emptiness
      plane = init_body_empty(4, 2, 0.0, 0.6, 0.4);
      create_cube!(plane, 20.0, 20.0, 0.2)
      plane.sv.x[3] = -3.0;

      window = init_screen();
      origo = origin();
      bodyvis = cube_vis(body);
      planevis = cube_vis(plane);
      pop!(planevis); # Remove local CS from plane

      gl._view(origo, window, camera=:perspective);
      gl._view(bodyvis, window, camera=:perspective);
      gl._view(planevis, window, camera=:perspective);

      sleep(1)
      t = 0.0;
      while t < simTime
            clear_forces!(body);
            vecadd!(body.f.F, g, body.md.m)
            body_dynamics!(body, delta_t)

            gl.set_arg!(bodyvis, :model, transformation(body))
            t += delta_t;
      end
      return nothing
end
