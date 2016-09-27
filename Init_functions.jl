# Functions used to initialize the simulation
#######################################
# Laskee sylinterin R_sektorin & Θ_sektorin.
function func_sektorit_cylinder!(R_sektori, Θ_sektori, n::Int64, k::Int64, Radius::Float64)
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
      local wi::Float64 = w / (k-1)

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
# Laskee kuution lokaalit karteesiset koordinaatit
function func_XYZ_cube_rigidbody!(XYZ, h, w, l, n, k)
      local hi::Float64 = h / 2.0
      local wi::Float64 = w / 2.0
      local li::Float64 = l / (k-1.0)

      for j in range(1,k)
            @inbounds begin
            XYZ[1,j][1] = wi
            XYZ[2,j][1] = wi
            XYZ[3,j][1] = -wi
            XYZ[4,j][1] = -wi
            XYZ[5,j][1] = wi

            XYZ[1,j][2] = -hi
            XYZ[2,j][2] = hi
            XYZ[3,j][2] = hi
            XYZ[4,j][2] = -hi
            XYZ[5,j][2] = -hi

            if j == 1
                  @simd for i in range(1,n+1)
                        XYZ[i,j][3] = -l/2.0
                  end
            else
                  @simd for i in range(1,n+1)
                        XYZ[i,j][3] = XYZ[i,1][3] + li*(j-1.0)
                  end
            end
            end #inbounds end
      end

      return nothing
end
