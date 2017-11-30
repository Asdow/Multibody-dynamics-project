# Define own mesh primitives/shapes here. Not using GeometryTypes primitives because they contain too many faces.
include("Cube.jl")
include("Sphere.jl")
include("Plane.jl")
#--------------------
"""
    init_shape(x::T, y::T, z::T) where {T<:Real}
Creates a cuboid Shape.
"""
function init_shape(x::T, y::T, z::T) where {T<:Real}
      body = cube_mesh(x, y, z)
      coll = cube_coll(x, y, z)
      bb = AABB(mpoint3D(zeros(T,3)), mpoint3D(zeros(T,3)))
      sh = Shape(body, coll, bb, deepcopy(bb) )
end
"""
    init_shape(radius::T) where {T<:Real}
Creates a sphere Shape.
"""
function init_shape(radius::T, center = point3D{T}(T(0.0),T(0.0),T(0.0))) where {T<:Real}
      body = sphere_mesh(radius)
      coll = CollShape(CollPrimitive{T}[])
      addSphere!(coll, radius, center)
      # coll = sphere_coll(radius, center)
      bb = AABB(mpoint3D(zeros(T,3)), mpoint3D(zeros(T,3)))
      sh = Shape(body, coll, bb, deepcopy(bb) )
end
