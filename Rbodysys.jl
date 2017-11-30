mutable struct RBodySystem{T<:Real}
    bodies::Array{RigidBody{T},1} # A list of bodies
    nb::Int64 # Number of bodies
    t::Float64 # time
end
function init_RSys(bodies::Array{RigidBody{T},1}) where {T}
    Rsys = RBodySystem(bodies, length(bodies), 0.0)
end
#-------------------------------------------------------
# Interface functions
getbodies(a::RBodySystem) = a.bodies
getnb(a::RBodySystem) = a.nb

import Base.length
length(a::RBodySystem) = a.nb
eltype(a::RBodySystem{T}) where {T} = T
eltype(a::Type{RBodySystem{T}}) where {T} = T

randomize!(a::RBodySystem) = randomize!(getbodies(a), length(a))
