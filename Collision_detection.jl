# Functions  for collision detection go here
struct Contact{T<:Real}
      A::Kappale # Kappale joka sisältää verteksin
      B::Kappale # Kappale joka sisältää pinnan

      p::sa.SVector{3,T} # Kontaktipisteen koordinaatit globaalissa koordinaatistossa
      n::sa.SVector{3,T} # Pinnan normaali, joka osoittaa ulospäin
      eA::sa.SVector{3,T} # Reunan suuntavektori kappaleelle A
      eB::sa.SVector{3,T} # Reunan suuntavektori kappaleelle B

      vf::Bool # true jos kontakti on verteksi/pinta
end
"""
    AABBb!(body::Kappale)
Goes through body's local coords and updates its local AABB.
"""
function AABBb!(body::Kappale)
      AABB!(body.sh.AABBb, body.sh.mesh.vertices)
      return nothing
end
"""
    AABB!(body::Kappale)
Goes through body's global coords and updates its global AABB.
"""
function AABB!(body::Kappale)
      xyz = [mpoint3D{Float64}(zeros(Float64,3)) for i in 1:size(body.sh.mesh.vertices,1)];
      body2world!(xyz, body.sv.x, body.aux.R, body.sh.mesh.vertices)
      AABB!(body.sh.AABB, xyz)
      return nothing
end
"""
    AABB!(bb::AABB, XYZ::Array{gt.Point3{T}}) where T
Goes through array's points and calculates the AABB inplace.
"""
function AABB!(bb::AABB, XYZ)
      @inbounds begin
            # Set first vertex' coords as AABB
            for i in 1:3
                  bb.min[i] = eltype(bb.max)(XYZ[1][i])
                  bb.max[i] = eltype(bb.max)(XYZ[1][i])
            end
            for i in eachindex(XYZ)
                  for j in 1:3
                        if XYZ[i][j] > bb.max[j]
                              bb.max[j] = eltype(bb.max)(XYZ[i][j])
                        end
                        if XYZ[i][j] < bb.min[j]
                              bb.min[j] = eltype(bb.max)(XYZ[i][j])
                        end
                  end
            end
      end
      return nothing
end
"""
    AABBOverlap(a::Kappale, b::Kappale)
Returns true if two bodies' AABBs overlap.
"""
function AABBOverlap(a::Kappale, b::Kappale)
      AABBOverlap(a.sh.AABB, b.sh.AABB)
end
"""
    AABBOverlap(a::AABB, b::AABB)
Returns true if two AABBs overlap.
"""
function AABBOverlap(a::AABB, b::AABB)
      d1 = mpoint3D(zeros(eltype(a.min),3))
      d2 = mpoint3D(zeros(eltype(a.min),3))
      vecminus!(d1, a.min, b.max)
      vecminus!(d2, b.min, a.max)
      if (d1.x > 0.0) || (d1.y > 0.0) || (d1.z > 0.0)
            return false
      elseif (d2.x > 0.0) || (d2.y > 0.0) || (d2.z > 0.0)
            return false
      end
      return true
end
