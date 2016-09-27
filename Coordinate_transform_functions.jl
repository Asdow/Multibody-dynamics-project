# Functions related to coordinate transformations
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
