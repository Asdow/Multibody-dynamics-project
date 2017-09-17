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
function MoI_cube!(J_array, mass, xl, yl, zl)
      fill!(J_array, 0.0)
      @inbounds begin
            J_array[1,1] = (mass*(zl^2 + yl^2)) / 12.0
            J_array[2,2] = (mass*(zl^2 + xl^2)) / 12.0
            J_array[3,3] = (mass*(yl^2 + xl^2)) / 12.0
      end
      return nothing
end
