# File for testing and interactive development
#######################################
include("Multiprocess.jl")
include("Packages.jl")
include("Shapes.jl")
include("Globals.jl")
include("apufunktiot.jl")
include("Quaternions.jl")
include("Datatypes.jl")
include("Particles.jl")
include("Visualization.jl")
include("Body_Init_functions.jl")
include("Moment_of_inertia_functions.jl")
include("Solids_functions.jl")
include("Collision_detection.jl")
include("Coordinate_transformations.jl")
include("Dynamics_functions.jl")
include("EOMs.jl")

function test()
      Rsys = init_RSys([create_cube(ones(3)..., 15.0, 0.6, 0.45) for i in 1:5]);
      nb = length(Rsys);
      bodies = Rsys.bodies;
      # Create plane so body doesn't fall into emptiness
      plane = create_cube(10.0,10.0,0.2, 0.0, 0.6, 0.45);
      plane.sv.x[3] = -3.0;

      window = init_screen();
      origo = origin();
      bodyvis = [body_vis(bodies[i]) for i in 1:nb];
      planevis = body_vis(plane);
      pop!(planevis); # Remove local CS from plane

      gl._view(origo, window, camera=:perspective);
      gl._view(bodyvis, window, camera=:perspective);
      gl._view(planevis, window, camera=:perspective);

      sleep(1)
      t = 0.0;
      while t < simTime
            clear_forces!(body);
            vecadd!(body.f.F, g, body.md.m)
            body_dynamics!(body, Δt)

            gl.set_arg!(bodyvis, :model, transformation(body))
            t += Δt;
            sleep(0.001);
      end
      return nothing
end

function testparticle()
      ps = ParticleSys([init_particle(rand()) for i in 1:10], 10)
      for i in 1:ps.n
            ps.p[i].x[:] = rand(3)
      end
      window = init_screen();
      origo = origin();

      return nothing
end
