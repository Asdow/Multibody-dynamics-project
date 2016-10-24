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

# Transforms body's local coordinates to world coordinates
function func_body_to_world(kpl::Kappale)
      @simd for j in range(1,kpl.pd.k)
          @simd for i in range(1,kpl.pd.n+1)
             @inbounds begin
                  kpl.world.XYZ[i,j] = kpl.sv.x + kpl.aux.Rot_mat*kpl.body.XYZ[i,j]
                  kpl.world.X[i,j] = kpl.world.XYZ[i,j][1]
                  kpl.world.Y[i,j] = kpl.world.XYZ[i,j][2]
                  kpl.world.Z[i,j] = kpl.world.XYZ[i,j][3]
             end
          end
      end
      return nothing
end

# Transforms body's local inertia to world
function func_inertia_body_to_world!(kpl::Kappale)
      #Inertia_world = Rotation_matrix_body * Inertia_body * Rotation_matrix_body'
      #a = @MMatrix(zeros(3,3))
      a = zeros(3,3)
      # body to world
      A_mul_B!(a, kpl.aux.Rot_mat, kpl.J.body)
      A_mul_B!(kpl.J.world, a, kpl.aux.Rot_mat')
      # body_inv to world_inv
      A_mul_B!(a, kpl.aux.Rot_mat, kpl.J.body_inv)
      A_mul_B!(kpl.J.world_inv, a, kpl.aux.Rot_mat')
      return nothing
end
