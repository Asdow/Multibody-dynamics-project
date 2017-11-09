# Functions used to calculate mass moment of inertias a
#######################################
"""
    MoI_cylinder!(J, mass, radius, w)
Calculates the moment of inertia of a cylinder. Updates it to J.
"""
function MoI_cylinder!(J, mass, radius, w)
    @assert size(J) == (3,3)
    fill!(J, 0.0)
    @inbounds begin
        J[1,1] = (mass/12.0) * (3.0*radius^2 + w^2)
        J[2,2] = copy(J[1,1])
        J[3,3] = (mass*radius^2) / 2.0
    end
    return nothing
end
"""
    MoI_cube!(J, mass, xl, yl, zl)
Calculates the moment of inertia of a cuboid, with sides xl, yl, & zl. Updates it to J.
"""
function MoI_cube!(J, mass, xl, yl, zl)
    @assert size(J) == (3,3)
    fill!(J, 0.0)
      @inbounds begin
            J[1,1] = (mass*(zl^2 + yl^2)) / 12.0
            J[2,2] = (mass*(zl^2 + xl^2)) / 12.0
            J[3,3] = (mass*(yl^2 + xl^2)) / 12.0
      end
      return nothing
end
"""
    MoI_sphere!(J, mass, radius)
Calculates the moment of inertia of a solid sphere. Updates it to J.
"""
function MoI_sphere!(J, mass, radius)
    @assert size(J) == (3,3)
    fill!(J, 0.0)
    @inbounds begin
        for i in 1:3
            J[i,i] = 2.0/5.0 * mass*radius^2
        end
    end
    return nothing
end
