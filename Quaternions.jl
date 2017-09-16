# Functions related to quaternions / Euler parameters
#######################################
# Mutable quaternions
mutable struct mQuaternion{T<:Real} <: sa.FieldVector{4, T}
      s::T
      v1::T
      v2::T
      v3::T
end
mQuaternion(a::Vector) = mQuaternion(a[1], a[2], a[3], a[4])
"""
    rotmat!(R::StaticArrays.MArray{Tuple{3,3},Float64,2,9}, q::mQuaternion)
Calculates the rotation matrix based on quaternions.
"""
function rotmat!(R::StaticArrays.MArray{Tuple{3,3},Float64,2,9}, q::mQuaternion)
      @inbounds begin
            e1_toiseen = q.s*q.s
            e2_toiseen = q.v1*q.v1
            e3_toiseen = q.v2*q.v2
            e4_toiseen = q.v3*q.v3
            e2e3 = q.v1*q.v2
            e1e4 = q.s*q.v3
            e2e4 = q.v1*q.v3
            e1e3 = q.s*q.v2
            e3e4 = q.v2*q.v3
            e1e2 = q.s*q.v1
            # Wikipedia
            R[1,1]=e1_toiseen + e2_toiseen - e3_toiseen - e4_toiseen
            R[2,1]=2.0*(e2e3 + e1e4)
            R[3,1]=2.0*(e2e4 - e1e3)

            R[1,2]=2.0*(e2e3 - e1e4)
            R[2,2]=e1_toiseen - e2_toiseen + e3_toiseen - e4_toiseen
            R[3,2]=2.0*(e3e4 + e1e2)

            R[1,3]=2.0*(e2e4 + e1e3)
            R[2,3]=2.0*(e3e4 - e1e2)
            R[3,3]=e1_toiseen - e2_toiseen - e3_toiseen + e4_toiseen
      end
      return nothing
end

"""
    q_vel!(q̇::mQuaternion, q::mQuaternion, ωs)
Calculates q̇ inplace.
"""
function q_vel!(q̇::mQuaternion, q::mQuaternion, ωs)
      for i in 1:4
            @inbounds q̇[i] = 0.5*(ωs[i,1]*q[1] + ωs[i,2]*q[2] + ωs[i,3]*q[3] + ωs[i,4]*q[4])
      end
      return nothing
end
"""
    q!(q::mQuaternion, q̇::mQuaternion, delta_t::Real)
Calculates q inplace.
"""
function q!(q::mQuaternion, q̇::mQuaternion, delta_t::Real)
      for i in 1:4
            @inbounds q[i] += q̇[i]*delta_t
      end
      return nothing
end


# Base functions
import Base: show, abs, normalize!, *

function show(io::IO, q::mQuaternion)
      print(io, "scalar: ",q.s, " vector: ",q.v1, "i ", q.v2, "j ", q.v3, "k")
end
abs(q::mQuaternion) = sqrt(q.s*q.s + q.v1*q.v1 + q.v2*q.v2 + q.v3*q.v3)
function normalize!(q::mQuaternion)
    pituus = abs(q)
    for i in 1:4
        @inbounds q[i] = q[i]/pituus
    end
    return nothing
end

(*)(q::mQuaternion, w::mQuaternion) = mQuaternion(q.s*w.s - q.v1*w.v1 - q.v2*w.v2 - q.v3*w.v3,
                                               q.s*w.v1 + q.v1*w.s + q.v2*w.v3 - q.v3*w.v2,
                                               q.s*w.v2 - q.v1*w.v3 + q.v2*w.s + q.v3*w.v1,
                                               q.s*w.v3 + q.v1*w.v2 - q.v2*w.v1 + q.v3*w.s)
