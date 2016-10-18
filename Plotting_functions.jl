# Functions related to plotting the simulation s
#######################################
using PyCall
using PyPlot
# Function to initialize the plot
function func_init_plot(xlim, ylim, zlim)
      fig = plt[:figure]()
      ax = fig[:add_subplot](111, projection="3d")
      # Määritellään koordinaatiston rajat
      ax[:set_xlim3d](-xlim,xlim)
      ax[:set_zlim3d](-zlim,zlim)
      ax[:set_ylim3d](-ylim,ylim)

      ax[:view_init](elev=92, azim=-90)
      fig[:canvas][:draw]()
      return fig, ax
end

# Function to plot the origin axes
function func_origin_axes(ax)
      X_akseli = [0.4 0.0 0.0; 0.0 0.0 0.0]
      Y_akseli = [0.0 0.0 0.0; 0.0 0.4 0.0]
      Z_akseli = [0.0 0.0 0.0; 0.0 0.0 0.4]
      origin_axes = Array(PyCall.PyObject, 3)
      origin_axes[1] = ax[:plot_wireframe](X_akseli[:,1], X_akseli[:,2], X_akseli[:,3], color="r")
      origin_axes[2] = ax[:plot_wireframe](Y_akseli[:,1], Y_akseli[:,2], Y_akseli[:,3], color="g")
      origin_axes[3] = ax[:plot_wireframe](Z_akseli[:,1], Z_akseli[:,2], Z_akseli[:,3], color="b")
      return origin_axes
end
function func_origin_axes!(ax, origin_axes)
      X_akseli = [0.4 0.0 0.0; 0.0 0.0 0.0]
      Y_akseli = [0.0 0.0 0.0; 0.0 0.4 0.0]
      Z_akseli = [0.0 0.0 0.0; 0.0 0.0 0.4]
      #origin_axes = Array(PyCall.PyObject, 3)
      origin_axes[1] = ax[:plot_wireframe](X_akseli[:,1], X_akseli[:,2], X_akseli[:,3], color="r")
      origin_axes[2] = ax[:plot_wireframe](Y_akseli[:,1], Y_akseli[:,2], Y_akseli[:,3], color="g")
      origin_axes[3] = ax[:plot_wireframe](Z_akseli[:,1], Z_akseli[:,2], Z_akseli[:,3], color="b")
      return nothing
end

# Function to plot the body's local axes
function func_origin_axes(ax, body)
      X_akseli = body.sv.x + body.aux.Rot_mat*[0.2; 0.0; 0.0]
      Y_akseli = body.sv.x + body.aux.Rot_mat*[0.0; 0.2; 0.0]
      Z_akseli = body.sv.x + body.aux.Rot_mat*[0.0; 0.0; 0.2]

      body_axes = Array(PyCall.PyObject, 3)
      body_axes[1] = ax[:plot_wireframe]( [X_akseli[1]; body.sv.x[1]], [X_akseli[2]; body.sv.x[2]], [X_akseli[3]; body.sv.x[3]], color="r" )
      body_axes[2] = ax[:plot_wireframe]([Y_akseli[1]; body.sv.x[1]], [Y_akseli[2]; body.sv.x[2]], [Y_akseli[3]; body.sv.x[3]], color="g")
      body_axes[3] = ax[:plot_wireframe]([Z_akseli[1]; body.sv.x[1]], [Z_akseli[2]; body.sv.x[2]], [Z_akseli[3]; body.sv.x[3]], color="b")
      return body_axes
end
function func_origin_axes!(ax, body, body_axes)
      X_akseli = body.sv.x + body.aux.Rot_mat*[0.2; 0.0; 0.0]
      Y_akseli = body.sv.x + body.aux.Rot_mat*[0.0; 0.2; 0.0]
      Z_akseli = body.sv.x + body.aux.Rot_mat*[0.0; 0.0; 0.2]

      #body_axes = Array(PyCall.PyObject, 3)
      body_axes[1] = ax[:plot_wireframe]( [X_akseli[1]; body.sv.x[1]], [X_akseli[2]; body.sv.x[2]], [X_akseli[3]; body.sv.x[3]], color="r" )
      body_axes[2] = ax[:plot_wireframe]([Y_akseli[1]; body.sv.x[1]], [Y_akseli[2]; body.sv.x[2]], [Y_akseli[3]; body.sv.x[3]], color="g")
      body_axes[3] = ax[:plot_wireframe]([Z_akseli[1]; body.sv.x[1]], [Z_akseli[2]; body.sv.x[2]], [Z_akseli[3]; body.sv.x[3]], color="b")
      return nothing
end

# Function to plot all the bodies
function func_draw_bodies!(ax, body_plot, Nbodies)
      for i in range(1,length(Nbodies))
            body_plot[i] = ax[:plot_wireframe](Nbodies[i].world.X, Nbodies[i].world.Y, Nbodies[i].world.Z)
      end
      return nothing
end
