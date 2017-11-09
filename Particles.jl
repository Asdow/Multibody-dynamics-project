# Particle system
abstract type Partikkeli end
struct Particle{T<:Real} <: Partikkeli
      x::mpoint3D{T} # Particle location in global coordinates
      ẋ::mpoint3D{T} # Particle velocity in global coordinates
      F::sa.MVector{3,T} # Forces effecting the particle in global coordinates
      m::T # Mass of the particle
end
struct ParticleSys{T<:Real}
      p::Vector{Particle{T}}; # array of particles
      n::Int64; # number of particles
      function ParticleSys(p::Vector{Particle{T}}) where {T}
            psys = ParticleSys(p, length(p))
      end
end
import Base: length;
length(p::ParticleSys) = p.n
#--------------------------------------
"""
    init_particle(m::T) where {T<:Real}
Creates a particle with mass m.
"""
function init_particle(m::T) where {T<:Real}
      p = Particle(mpoint3D(zeros(T,3)), mpoint3D(zeros(T,3)), sa.MVector{3, T}(zeros(T,3)), m)
end
#particlemsh = gl.GLNormalMesh(gl.loadasset("I:/Julia/sphere_radius0.1.obj"))
