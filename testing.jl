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
func_body_to_world(cube)
# Simulation times
t_init = 0.0
delta_t = 0.005
t_end = 0.50
cube.torque_body[1] = 0.0
cube.torque_body[2] = 0.0
cube.torque_body[3] = 0.0

for t = collect(t_init:delta_t:t_end)
      cube.torque_world = cube.Rot_mat*cube.torque_body

      cube.L = cube.L + cube.torque_world*delta_t
      cube.omega_world = cube.J_world_inv*cube.L
      func_omega_tilde4!(cube.omega_tilde4, cube.omega_world)
      func_euler_vel!(cube.qdot, cube.q, cube.omega_tilde4)
      func_euler!(cube.qa, cube.qdot, cube.q, delta_t)
      func_euler_normalize!(cube.q)
      func_rotation_matrix!(cube.Rot_mat,cube.q)
      cube.J_world_inv = cube.Rot_mat*cube.J_body_inv*(cube.Rot_mat')

      # Update world coordinates
      func_body_to_world(cube)
end
using GLVisualize, GeometryTypes, Reactive, GLAbstraction, Colors,GLWindow, ModernGL
"""
This code should be executed only one time per julia session!!
If you accidantly close the window, you can call this again.
"""
function init_GLscreen(res=(800,600))
      # giving the window a transparent background color makes it transparent to
        # the previous frame. It's arguable, if that's really how things should be,
        # but that's how it currently works ;)
        window = glscreen(resolution=res, background=RGBA(0,0,0,0))
        timesignal = Signal(0.0)
        speed = Signal(1/60)

        #glClear(GL_COLOR_BUFFER_BIT)
        # this is equivalent to @async renderloop(window).
        # but now you can do other stuff before an image is rendered
        # the @async is used to make this non blocking for working in the REPL/Atom
        @async while isopen(window)
              # timesignal update value (here 0.01) doesn't matter for particle simulation, but controls the color update speed
             push!(timesignal, value(timesignal)+delta_t)
             render_frame(window)
             sleep(value(speed))
        end
        window, timesignal, speed
 end

window, timesignal, speed = init_GLscreen()
surface = visualize((cube.X_world,cube.Y_world,cube.Z_world), :surface)
_view(surface, window)
renderloop(window)
