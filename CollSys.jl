abstract type CollSys{T<:Real} end
abstract type Collshape{T<:Real} end
abstract type CollPrimitive{T<:Real} end

# CollPrimitive defines a collision primitive and its location in a body's local frame.
struct CollSphere{T} <: CollPrimitive{T}
    # Region R = { {x,y,z} | (x-c.x)^2 + (y-c.y)^2 + (z-c.z)^2 <= r^2 }
    c::point3D{T} # Sphere center. Defined in body local coordinates. If c = zeros(3), it coincides with the reference frame location.
    r::T # Radius of the sphere.
end

# CollShape holds a list of coll. primitives. A CollShape can be used to add several collision primitives to a body.
struct CollShape{T} <: Collshape{T}
    primitives::Array{CollPrimitive{T},1}
end
getprimitives(cs::Collshape{T}) where {T} = cs.primitives

"""
    addSphere!(cs::Collshape{T}, r::T, c::point3D{T}=point3D(zero(T),zero(T),zero(T))) where {T}
Adds a CollSphere to cs. If the sphere's center is not given, it defaults to local origin.
"""
function addSphere!(cs::Collshape{T}, r::T, c::point3D{T}=point3D(zero(T),zero(T),zero(T))) where {T}
    sphere = CollSphere(c,r)
    list = getprimitives(cs)
    push!(list, sphere)
    return nothing
end
# CollSystem holds a list of all the currently available CollShapes.
# TODO think of a way to remove CollShapes when no body is using them anymore.
struct CollSystem{T} <: CollSys{T}
    Collshapes::Array{CollShape{T},1}
end
function createCollSystem()
    Csys = CollSystem(CollShape{Float64}[])
end
getCollShapes(Csys::CollSys{T}) where {T} = Csys.Collshapes

"""
    addCollShape!(Csys::CollSys{T}, cshape = CollShape(CollPrimitive{T}[])) where {T}
Adds a CollShape to a Csys. If no cshape is defined, it adds an empty Collshape to Csys.
"""
function addCollShape!(Csys::CollSys{T}, cshape = CollShape(CollPrimitive{T}[])) where {T}
    # cshape = CollShape(CollPrimitive{T}[])
    list = getCollShapes(Csys)
    push!(list, cshape)
    return nothing
end
