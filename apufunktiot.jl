# Helper functions
"Funktio joka laskee 3D vektorien a ja b ristitulon ja tallentaa tuloksen vektoriin y"
function cross!(y, a, b)
      @inbounds begin
            y[1] = a[2]*b[3] - a[3]*b[2]
            y[2] = a[3]*b[1] - a[1]*b[3]
            y[3] = a[1]*b[2] - a[2]*b[1]
      end
      return nothing
end

"Funktio joka laskee vektorivähennyksen a - b ja tallentaa sen vektoriin y"
function vecminus!(y, a, b)
      l = length(y)
      @simd for i in 1:l
            @inbounds y[i] = a[i] - b[i]
      end
      return nothing
end

"""
    vecplus!(y, a, b)
Funktio joka laskee vektorisumman a + b ja tallentaa sen vektoriin y. Vektorien pituuksien oltava samat.
"""
function vecplus!(y, a, b)
      l = length(y)
      @simd for i in 1:l
            @inbounds y[i] = a[i] + b[i]
      end
      return nothing
end
"""
    vecadd!(a, b, c::Real)
Laskee a + b*c ja tallentaa sen vektoriin a. Vektorien pituuksien pitää olla 3.
"""
function vecadd!(a, b, c::Real)
      for i in 1:3
            @inbounds a[i] += b[i]*c
      end
      return nothing
end
"""
    vecadd2!(a, b, c::Real)
Laskee b*c ja tallentaa sen vektoriin a. Vektorien pituuksien pitää olla 3.
"""
function vecadd2!(a, b, c::Real)
      for i in 1:3
            @inbounds a[i] = b[i]*c
      end
      return nothing
end
"""
    inverse!(minv, m)
Calculates 3x3 matrix inverse inplace and saves it to minv.
"""
function inverse!(minv, m)
      @inbounds begin
            if iszero(m) == true
                  fill!(minv, zero(eltype(minv)))
            else
                  m2233_3223 = m[2,2] * m[3,3] - m[3,2] * m[2,3]
                  m2133 = m[2,1] * m[3,3]
                  m2331 = m[2,3] * m[3,1]
                  m2132_2231 = m[2,1] * m[3,2] -  m[2,2] * m[3,1]
                  # computes the inverse of a matrix m
                  invdet = 1.0 / (m[1,1] * m2233_3223 - m[1,2] * (m2133 - m2331) + m[1,3] * m2132_2231);
                  # inverse of matrix m
                  minv[1,1] = m2233_3223 * invdet;
                  minv[2,1] = (m2331 - m2133) * invdet;
                  minv[3,1] = m2132_2231 * invdet;
                  minv[1,2] = (m[1,3] * m[3,2] - m[1,2] * m[3,3]) * invdet;
                  minv[2,2] = (m[1,1] * m[3,3] - m[1,3] * m[3,1]) * invdet;
                  minv[3,2] = (m[3,1] * m[1,2] - m[1,1] * m[3,2]) * invdet;
                  minv[1,3] = (m[1,2] * m[2,3] - m[1,3] * m[2,2]) * invdet;
                  minv[2,3] = (m[2,1] * m[1,3] - m[1,1] * m[2,3]) * invdet;
                  minv[3,3] = (m[1,1] * m[2,2] - m[2,1] * m[1,2]) * invdet;
            end
      end
      return nothing
end
