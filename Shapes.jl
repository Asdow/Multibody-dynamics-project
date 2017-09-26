# Define own mesh primitives/shapes here. Not using GeometryTypes primitives because they contain too many faces.
"""
    cube(x::T, y::T, z::T) where {T<:Integer}
Creates a GLNormalMesh cube primitive with sides x, y and z.
"""
# NOTE Horrible hack. Need a better way of doing this.
function cube(x::T, y::T, z::T) where {T<:Real}
      vertices = gt.Point{3,Float32}[]
      w = x/2.0
      h = y/2.0
      l = z/2.0
      push!(vertices, gt.Point3f0(-w, -h, l))
      push!(vertices, gt.Point3f0(-w, h, l))
      push!(vertices, gt.Point3f0(-w, -h, -l))
      push!(vertices, gt.Point3f0(-w, h, -l))
      push!(vertices, gt.Point3f0(w, -h, l))
      push!(vertices, gt.Point3f0(w, h, l))
      push!(vertices, gt.Point3f0(w, -h, -l))
      push!(vertices, gt.Point3f0(w, h, -l))
      faces = gt.Face{3, gt.OffsetInteger{-1,UInt32}}[]
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [1,2,4]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [1,4,3]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [3,4,8]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [3,8,7]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [7,8,6]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [7,6,5]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [5,6,2]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [5,2,1]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [3,7,5]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [3,5,1]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [8,4,2]])
      push!(faces, [gt.OffsetInteger{-1,UInt32}(i) for i in [8,2,6]])
      normals = GeometryTypes.Normal{3,Float32}[]
      x1 = 0.816497;
      x2 = 0.333333;
      y1 = 0.408248;
      y2 = 0.666667;
      z1 = y1;
      z2 = y2;
      push!(normals, gt.Normal{3,Float32}(-x1, -y1, z1))
      push!(normals, gt.Normal{3,Float32}(-x2, y2, z2))
      push!(normals, gt.Normal{3,Float32}(-x2, -y2, -z2))
      push!(normals, gt.Normal{3,Float32}(-x1, y1, -z1))
      push!(normals, gt.Normal{3,Float32}(x2, -y2, z2))
      push!(normals, gt.Normal{3,Float32}(x1, y1, z1))
      push!(normals, gt.Normal{3,Float32}(x1, -y1, -z1))
      push!(normals, gt.Normal{3,Float32}(x2, y2, -z2))
      msh = gl.GLNormalMesh(vertices, faces, normals, Void[], Void(), Void(), Void[])
      return msh
end
