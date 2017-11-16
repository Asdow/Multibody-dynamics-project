# functions and structs related to spheres.#
############################################
"""
    sphere_mesh(r)
Creates a GLNormalMesh sphere primitive with radius r.
"""
function sphere_mesh(r)
    sphere = gl.Sphere{Float32}(gl.Point3f0(0), Float32(r));
    vertices = gl.decompose(gl.Point3f0, sphere, 18);
    faces = gl.decompose(gl.GLTriangle, sphere, 18);
    msh = gl.GLNormalMesh(vertices=vertices, faces=faces)
end
"""
    sphere_coll(radius::T, center = point3D{T}(T(0.0),T(0.0),T(0.0))) where {T<:Real}
Create a sphere collision shape.
"""
function sphere_coll(radius::T, center = point3D{T}(T(0.0),T(0.0),T(0.0))) where {T<:Real}
    CollSphere(center,radius)
end
"""
    getsphereworldpos(s::CollSphere, b::Kappale)
Get the world position of a CollSphere.
"""
function getsphereworldpos(s::CollSphere, b::Kappale)
    center = b.sv.x + b.aux.R*s.c # Transform sphere center to world.
end
getsphereworldpos(b::Kappale) = getsphereworldpos(b.sh.coll, b)
