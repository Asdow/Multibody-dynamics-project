# Functions used to calculate mass moment of inertias a
#######################################
# Calculates the moment of inertia of a cylinder
function func_MoI_cylinder!(J_array, mass, radius, w)
      J_array[1,1] = (mass/12.0) * (3.0*radius^2 + w^2)
      J_array[2,2] = copy(J_array[1,1])
      J_array[3,3] = (mass*radius^2) / 2.0

      return nothing
end
# Calculates the moment of inertia of a cube
function func_MoI_cube!(J_array, mass, h, w, l)
      J_array[1,1] = (mass*(h^2 + l^2)) / 12.0
      J_array[2,2] = (mass*(l^2 + w^2)) / 12.0
      J_array[3,3] = (mass*(w^2 + h^2)) / 12.0

      return nothing
end
