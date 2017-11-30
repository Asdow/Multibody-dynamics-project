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


"""
    TestSphereHalfspace(s::CollSphere, p::Plane, b::Kappale)
Determine whether sphere s intersects the negative halfspace of plane p. Body is needed to transform sphere center location from local coordinates to inertial frame.
"""
function TestSphereHalfspace(s::CollSphere, p::Plane, b::Kappale)
    center = getsphereworldpos(s, b) # Transform sphere center to world.
    dist = dot(center, p.normal) - p.d; # Check for intersection.
    return dist <= s.r
end
TestSphereHalfspace(p::Plane, b::Kappale) = TestSphereHalfspace(b.sh.coll, p, b)
TestSphereHalfspace(b::Kappale, p::Plane) = TestSphereHalfspace(p, b)

function getcollpos(s::CollSphere, p::Plane, b::Kappale)
    center = getsphereworldpos(s,b) # Transform sphere center to world.
    dist = dot(center, p.normal) - p.d; # Check for intersection.
    point = point3D(center - dist*p.normal) # Contact point is on the plane.
end
getcollpos(b::Kappale, p::Plane) = getcollpos(b.sh.coll, p, b)

function getcollposdepth(s::CollSphere, p::Plane, b::Kappale)
    center = getsphereworldpos(s,b) # Transform sphere center to world.
    dist = dot(center, p.normal) - p.d; # Check for intersection.
    point = point3D(center - dist*p.normal) # Contact point is on the plane.
    depth = s.r-dist
    return point, depth
end
getcollposdepth(b::Kappale, p::Plane) = getcollposdepth(b.sh.coll, p, b)

function getcollposdepth(as::CollSphere, bs::CollSphere, a::Kappale, b::Kappale)
    # Transform spheres to world frame
    aw = getsphereworldpos(as,a)
    bw = getsphereworldpos(bs,b)
    # Vector from B to A
    BA = aw-bw
    # unit vector pointing to sphere A
    n̂ = normalize(BA)

    depth = (as.r+bs.r) - norm(BA)
    # Coll. point is the average of coll. depth along coll. normal
    point = 0.5*( aw+as.r*n̂ + bw-bs.r*n̂ )
    return point, depth, n̂
end
getcollposdepth(a::Kappale, b::Kappale) = getcollposdepth(a.sh.coll, b.sh.coll, a, b)

function checkcoll(p::Plane, bodies::Array{T,1}, nb=length(bodies)) where {T<:Kappale}
    # Go through the list of bodies and check if any collide with plane. If true, add body index to list.
    list = Int64[];
    for i in 1:nb
        if TestSphereHalfspace(bodies[i], p)
            push!(list,i)
        end
    end
    return list
end
checkcoll(Rsys::RBodySystem, p::Plane) = checkcoll(p, Rsys.bodies, Rsys.nb)
checkcoll(p::Plane, Rsys::RBodySystem) = checkcoll(Rsys, p)

function checkcoll(bodies::Array{T,1}, nb=length(bodies)) where {T<:Kappale}
    # Go through the list of bodies and check if any collide with each other. If true, add body indices to list.
    list = Tuple{Int64,Int64}[];
    for j in 1:(nb-1)
        for i in (j+1):nb
            if TestSphereSphere(bodies[j], bodies[i])
                push!(list,(j,i))
            end
        end
    end
    return list
end
checkcoll(Rsys::RBodySystem) = checkcoll(Rsys.bodies, Rsys.nb)

function TestSphereSphere(a::CollSphere, b::CollSphere, ab::Kappale, bb::Kappale)
    # Transform spheres to world frame
    aw = getsphereworldpos(a,ab)
    bw = getsphereworldpos(b,bb)
    # Calculate squared distance between centers
    d = aw - bw;
    dist2 = dot(d, d);
    # Spheres intersect if squared distance is less than squared sum of radii
    radiusSum = a.r + b.r;
    return dist2 <= radiusSum * radiusSum;
end
TestSphereSphere(a::RigidBody{T}, b::RigidBody{T}) where {T} = TestSphereSphere(a.sh.coll, b.sh.coll, a, b)
