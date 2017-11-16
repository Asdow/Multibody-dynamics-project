# File for testing and interactive development
#######################################
include("Multiprocess.jl")
include("Packages.jl")
include("Globals.jl")
include("Quaternions.jl")
include("Datatypes.jl")
include("Shapes.jl")
include("apufunktiot.jl")
include("Visualization.jl")
include("Body_Init_functions.jl")
include("Moment_of_inertia_functions.jl")
include("Solids_functions.jl")
include("Collision_detection.jl")
include("Contact_modeling.jl")
include("Coordinate_transformations.jl")
include("Dynamics_functions.jl")
include("EOMs.jl")

function test()
    # Rsys = init_RSys([create_cube(ones(3)..., 15.0, 0.6, 0.45) for i in 1:5]);
    Rsys = init_RSys([create_sphere(1.0, 100.0, 0.6, 0.4) for i in 1:2]);
    nb = length(Rsys);
    bodies = Rsys.bodies;
    randomize!(Rsys)
    for i in 1:Rsys.nb
        Rsys.bodies[i].sv.x.z += 2.0;
    end
    # Create plane so body doesn't fall into emptiness
    plane = Plane([0.0,0.0, 1.0], point3D(0.0, 0.0, 0.0))

    window = init_screen();
    origo = origin();
    bodyvis = [body_vis(bodies[i]) for i in 1:nb];

    gl._view(origo, window, camera=:perspective);
    gl._view(bodyvis, window, camera=:perspective);

    sleep(1)
    while Rsys.t < simTime
        list = checkcoll(Rsys, plane);
        Flist = resolvecoll(Rsys, list, plane)

        eval_forces!(Rsys, Flist);
        dynamics!(Rsys)

        body_vis!(bodyvis, Rsys)
        Rsys.t += Δt;
        sleep(0.001);
    end
    return nothing
end
