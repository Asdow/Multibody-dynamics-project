# Functions used to initialize the different bodies
#######################################
"""
    create_cube(x::T1, y::T1, z::T1, m::T2, mus::T2, muk::T2) where {T1<:Real, T2<:Float64}
Create a cuboid shaped body. Sides are x,y,z, mass is m and friction coefficients are mus and muk.
"""
function create_cube(x::T1, y::T1, z::T1, m::T2, mus::T2, muk::T2) where {T1<:Real, T2<:Float64}
      sh = init_shape(x, y, z)
      md = init_MassData(m)
      sv = init_StateVariables()
      dv = init_Derivates()
      aux = init_Aux()
      f = init_Voimat()
      knt = init_PenMethod()
      kd = init_Lugre(mus, muk)
      body = RigidBody( sh, md, sv, dv, aux, f, knt, kd )
      AABBb!(body)
      rotmat!(body.aux.R, body.sv.q)
      inverse!(body.aux.R_inv, body.aux.R)
      # Moment of inertias
      MoI_cube!(body.md.Jb, body.md.m, x, y, z)
      inverse!(body.md.Jb_inv, body.md.Jb)
      Jb_inv2world!(body)
      Jb2world!(body)
      return body
end
"""
    cube_coords!(XYZ::Array{point3D{Float64},2}, h::T, w::T, l::T) where {T<:Float64}
Laskee kuution lokaalit karteesiset koordinaatit. Toimii vain jos Arrayn koko on 4x2
"""
function cube_coords!(XYZ::Array{point3D{Float64},2}, xl::T, yl::T, zl::T) where {T<:Float64}
      @assert size(XYZ) == (4,2) "XYZ size != (4,2)"
      n = size(XYZ,1)
      k = size(XYZ,2)
      local xi = xl / 2.0
      local yi = yl / 2.0
      local zi = zl / 2.0
      for j in 1:k
            @inbounds begin
                  XYZ[1,j].x = xi
                  XYZ[2,j].x = xi
                  XYZ[3,j].x = -xi
                  XYZ[4,j].x = -xi
                  XYZ[1,j].y = -yi
                  XYZ[2,j].y = yi
                  XYZ[3,j].y = yi
                  XYZ[4,j].y = -yi
            end
      end
      for i in 1:n
            XYZ[i,1].z = -zi
            XYZ[i,2].z = zi
      end
      return nothing
end
