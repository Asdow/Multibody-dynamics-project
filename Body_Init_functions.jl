# Functions used to initialize the different bodies
#######################################
# Laskee sylinterin R_sektorin & Θ_sektorin.
function func_sektorit_cylinder!(R_sektori::AbstractArray, Θ_sektori::AbstractArray, n::Integer, k::Integer, Radius::AbstractFloat)
    @simd for j in range(1,k)
        # Yhden sektorin kulma
        #Θ_sektori = 2*pi / n
        # Jaetaan profiili n:ään pisteeseen
        @simd for i in range(1,n+1)
            @inbounds R_sektori[i,j] = Radius# * (rand(89:105)/100)
            @inbounds Θ_sektori[i] = (2*pi / n) * (i-1)
        end
    end

    # Muokataan profiilin datan alkupisteitä, jotta plottauksessa profiili ei ole aukinainen
    @simd for i in range(1,k)
        @inbounds R_sektori[n+1,i] = copy(R_sektori[1,i])
    end

    return nothing
end
# Laskee sylinterin lokaalit karteesiset koordinaatit
function func_XYZ_cylinder_rigidbody!(XYZ, R_sektori, Θ_sektori, n, k, w)
      local wi = w / (k-1.0)

      @simd for j in range(1,k)
            @simd for i in 1:(n+1)
                  @inbounds XYZ[i,j][1],XYZ[i,j][2] = func_cartesis(R_sektori[i,j], Θ_sektori[i])
                  if j == 1
                        @inbounds XYZ[i,j][3] = -w/2.0
                  else
                        @inbounds XYZ[i,j][3] = XYZ[i,1][3] + wi*(j-1.0)
                  end
            end
      end
      return nothing
end
"""
    cube_rigidbody!(XYZ::Array{point3D{Float64},2}, h::T, w::T, l::T) where {T<:Float64}
Laskee kuution lokaalit karteesiset koordinaatit. Toimii vain jos Arrayn koko on 4x2
"""
function cube_rigidbody!(XYZ::Array{point3D{Float64},2}, h::T, w::T, l::T) where {T<:Float64}
      @assert size(XYZ) == (4,2) "XYZ size != (4,2)"
      n = size(XYZ,1)
      k = size(XYZ,2)
      local wi = w / 2.0
      local hi = h / 2.0
      local li = l / 2.0
      for j in 1:k
            @inbounds begin
                  XYZ[1,j].x = wi
                  XYZ[2,j].x = wi
                  XYZ[3,j].x = -wi
                  XYZ[4,j].x = -wi
                  XYZ[1,j].y = -hi
                  XYZ[2,j].y = hi
                  XYZ[3,j].y = hi
                  XYZ[4,j].y = -hi
            end
      end
      for i in 1:n
            XYZ[i,1].z = -li
            XYZ[i,2].z = li
      end
      return nothing
end
