struct Plane{T} <: Collshape{T}
    normal::sa.SVector{3,T} # Plane normal. Assumed to be a unit vector. Points x on the plane satisfy Dot(n,x) = d
    d::T # d = dot(n,p) for a given point p on the plane
end
function Plane(vector, point)
    n = sa.SVector{3,eltype(vector)}(normalize(vector))
    d = convert(eltype(vector), dot(n, point))
    Plane(n, d)
end
