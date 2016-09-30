# Functions to calculate volumes of bodies 
#######################################
# Volume of cylinder
function func_vol_cylinder(Radius::Float64, w::Float64)
      volume = pi*Radius^2*w
      return volume
end
# Volume of a cube
function func_vol_cube(h, w, l)
      volume = h*w*l
      return volume
end
