# Functions related to coordinate transformations s
#######################################
# Funktio joka muuttaa asteet radiaaneiksi
function func_radians(kulma::Float64)
    kulma * (pi/180.0)
end

# Funktio joka muuttaa radiaanit asteiksi
function func_asteet(kulma::Float64)
    kulma * (180.0/pi)
end

# Funktio joka muuttaa polaarikoordinaatin carteesiseen koordinaatistoon
function func_cartesis(R::Float64, Θ::Float64)
    x::Float64 = R * cos(Θ)
    y::Float64 = R * sin(Θ)
    return x,y
end

# Funktio joka muuttaa karteesiset koordinaatit polaariseen koordinaatistoon
function func_polar(X::Float64, Y::Float64)
    R = sqrt(X^2 + Y^2)
    Θ = atan(Y/X)
    return R,Θ
end

# Muutetaan kappaleen karteesiset koordinaatit polaariseen koordinaatistoon
function func_polar_array!(R::Array, Θ::Array, XYZ::Array, n, k)
      @simd for j in range(1,k)
            @simd for i in range(1,n+1)
                  @inbounds R[i,j], Θ[i] = func_polar(XYZ[i,j][1], XYZ[i,j][2])
            end
      end
    return nothing
end

# Muutetaan kappaleen karteesiset koordinaatit polaariseen koordinaatistoon
function func_cartesis_array!(R::Array, Θ::Array, XYZ::Array, n, k)
      @simd for j in range(1,k)
            @simd for i in range(1,n+1)
                  @inbounds XYZ[i,j][1], XYZ[i,j][2] = func_cartesis(R[i,j], Θ[i])
            end
      end
    return nothing
end

"""
    body2world!(kpl::Kappale)
Transforms body's local coordinates to world.
"""
function body2world!(kpl::Kappale)
      body2world!(kpl.coords.XYZ, kpl.sv.x, kpl.aux.R, kpl.coords.XYZb)
      return nothing
end
"""
    body2world!(XYZ, svx, R, XYZb)
Transforms local coords to world. Saves to XYZ.
"""
function body2world!(XYZ, svx, R, XYZb)
      for j in 1:size(XYZ,2)
           for i in 1:size(XYZ,1)
                @inbounds A_mul_B!(XYZ[i,j], R, XYZb[i,j])
                @inbounds XYZ[i,j] .+= svx
          end
      end
      return nothing
end
"""
    Jb_inv2world!(kpl::Kappale)
Transforms body's local inverted inertia matrix to world.
"""
function Jb_inv2world!(kpl::Kappale)
      J_local2global!(kpl.md.J_inv, kpl.aux.R, kpl.md.Jb_inv)
      return nothing
end
"""
    Jb2world!(kpl::Kappale)
Transforms body's local inertia matrix to world.
"""
function Jb2world!(kpl::Kappale)
      J_local2global!(kpl.md.J, kpl.aux.R, kpl.md.Jb)
      return nothing
end
"""
    J_local2global!(J, rotmat_body, Jbody)
Transforms 3x3 matrix from local to world. Saves to J.
"""
function J_local2global!(J, rotmat_body, Jbody)
      @inbounds begin
            RotJ11 = 0.0
            RotJ21 = 0.0
            RotJ31 = 0.0
            RotJ12 = 0.0
            RotJ22 = 0.0
            RotJ32 = 0.0
            RotJ13 = 0.0
            RotJ23 = 0.0
            RotJ33 = 0.0
            for i in 1:3
                  RotJ11 += rotmat_body[1,i]*Jbody[i,1]
                  RotJ21 += rotmat_body[2,i]*Jbody[i,1]
                  RotJ31 += rotmat_body[3,i]*Jbody[i,1]
                  RotJ12 += rotmat_body[1,i]*Jbody[i,2]
                  RotJ22 += rotmat_body[2,i]*Jbody[i,2]
                  RotJ32 += rotmat_body[3,i]*Jbody[i,2]
                  RotJ13 += rotmat_body[1,i]*Jbody[i,3]
                  RotJ23 += rotmat_body[2,i]*Jbody[i,3]
                  RotJ33 += rotmat_body[3,i]*Jbody[i,3]
            end
            J[1,1] = RotJ11*rotmat_body[1,1] + RotJ12*rotmat_body[1,2] + RotJ13*rotmat_body[1,3]
            J[2,1] = RotJ21*rotmat_body[1,1] + RotJ22*rotmat_body[1,2] + RotJ23*rotmat_body[1,3]
            J[3,1] = RotJ31*rotmat_body[1,1] + RotJ32*rotmat_body[1,2] + RotJ33*rotmat_body[1,3]

            J[1,2] = RotJ11*rotmat_body[2,1] + RotJ12*rotmat_body[2,2] + RotJ13*rotmat_body[2,3]
            J[2,2] = RotJ21*rotmat_body[2,1] + RotJ22*rotmat_body[2,2] + RotJ23*rotmat_body[2,3]
            J[3,2] = RotJ31*rotmat_body[2,1] + RotJ32*rotmat_body[2,2] + RotJ33*rotmat_body[2,3]

            J[1,3] = RotJ11*rotmat_body[3,1] + RotJ12*rotmat_body[3,2] + RotJ13*rotmat_body[3,3]
            J[2,3] = RotJ21*rotmat_body[3,1] + RotJ22*rotmat_body[3,2] + RotJ23*rotmat_body[3,3]
            J[3,3] = RotJ31*rotmat_body[3,1] + RotJ32*rotmat_body[3,2] + RotJ33*rotmat_body[3,3]
      end
      return nothing
end

"""
    rotmat!(body::Kappale)
Calculates body's rotation matrix from quaternions.
"""
function rotmat!(body::Kappale)
      rotmat!(body.aux.R, body.sv.q)
end
