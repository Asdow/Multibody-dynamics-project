# File for testing and interactive development
#######################################
workspace()
include("Datatypes.jl")
include("Init_functions.jl")
include("Solids_functions.jl")
include("Volume_functions.jl")
include("Moment_of_inertia_functions.jl")
include("Coordinate_transform_functions.jl")
include("Quaternion_functions.jl")


cube = func_init_body_empty(4, 2)
fill_rigidbody_cube!(cube, cube.n, cube.k, 1.0, 2.0, 3.0, 100.0)
func_polar_array!(cube.R, cube.Θ, cube.XYZ_body, cube.n, cube.k)
