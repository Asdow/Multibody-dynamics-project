# Functions related to plotting the simulation
#######################################
using Plots
glvisualize()
gr()

#############
n = 100
ts = linspace(0,8π,n)
x = ts .* map(cos,ts)
y = (0.1ts) .* map(sin,ts)
z = 1:n
plot(x,y,z,zcolor=reverse(z),m=(10,0.8,:blues,stroke(0)),leg=false,cbar=true,w=5)
plot!(zeros(n),zeros(n),1:n,w=10)
#####################
plot(Plots.fakedata(50,5),w=3)
################
n = 50
x = linspace(-3, 3, n)
y = x

z = Array(Float64, n, n)
f(x, y) = cos(x^2 + y^2) / (1 + x^2 + y^2)
for i in 1:n
    for j in 1:n
        z[j, i] = f(x[i], y[j])
    end
end

surface(x, y, z)

#########
x = Array(Float64, cube.n+1, cube.k)
y = Array(Float64,cube.n+1, cube.k)
z = Array(Float64,cube.n+1, cube.k)

cube.q[1] = -0.311
cube.q[2] = 0.544
cube.q[3] = 0.531
cube.q[4] = 0.571
func_rotation_matrix!(cube.Rot_mat,cube.q)

for j in range(1,cube.k)
      for i in range(1,cube.n+1)
            cube.XYZ_world[i,j] = cube.Rot_mat*cube.XYZ_body[i,j]
            x[i,j] = cube.XYZ_world[i,j][1]
            y[i,j] = cube.XYZ_world[i,j][2]
            z[i,j] = cube.XYZ_world[i,j][3]
      end
end
cube_graphics = surface(x,y,z);
#surface(vcat(x...), vcat(y...), vcat(z...))
gui()
