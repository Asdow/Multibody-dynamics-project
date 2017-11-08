# Functions for the rigid body dynamics
"""
    clear_forces!(kpl::T) where {T<:Kappale}
Clears body's forces.
"""
function clear_forces!(kpl::T) where {T<:Kappale}
    fill!(kpl.f.F, 0.0);
    fill!(kpl.f.Fb, 0.0);
    fill!(kpl.f.T, 0.0);
    fill!(kpl.f.Tb, 0.0);
    return nothing
end
function clear_forces!(bodies::Array{T,1}, nb=length(bodies)) where {T<:Kappale}
    for i in 1:nb
        clear_forces!(bodies[i])
    end
    return nothing
end
clear_forces!(Rsys::RBodySystem{T}) where {T} = clear_forces!(Rsys.bodies, Rsys.nb)

"""
    skew4!(ωskew, ω)
Calculates 4x4 skew symmetric matrix in global coordinate frame
"""
function skew4!(ωskew, omega)
    @assert size(ωskew) == (4,4)
    @assert size(omega) == (3,)
    @inbounds begin
          ωskew[1,1] = 0.0
          ωskew[2,1] = omega[1]
          ωskew[3,1] = omega[2]
          ωskew[4,1] = omega[3]

          ωskew[1,2] = -omega[1]
          ωskew[2,2] = 0.0
          ωskew[3,2] = omega[3]
          ωskew[4,2] = -omega[2]

          ωskew[1,3] = -omega[2]
          ωskew[2,3] = -omega[3]
          ωskew[3,3] = 0.0
          ωskew[4,3] = omega[1]

          ωskew[1,4] = -omega[3]
          ωskew[2,4] = omega[2]
          ωskew[3,4] = -omega[1]
          ωskew[4,4] = 0.0
    end
    return nothing
end
"""
    body_translation!(kpl::Kappale, delta_t::Real)
Calculates body's linear momentum, linear velocity and reference point global location. Updates all three inplace.
"""
function body_translation!(kpl::Kappale, delta_t::Real)
      @inbounds begin
            # Linear momentum
            vecadd!(kpl.sv.P, kpl.f.F, delta_t)
            # Massakeskipisteen nopeus
            vecadd2!(kpl.dv.ẋ, kpl.sv.P, kpl.md.m_inv)
            # Massakeskipisteen asema
            vecadd!(kpl.sv.x, kpl.dv.ẋ, delta_t)
      end
      return nothing
end
"""
    body_rotation!(kpl::Kappale, delta_t::Real)
Calculates body's angular momentum, angular velocity, reference point global orientation. Updates all inplace.
"""
function body_rotation!(kpl::Kappale, delta_t::Float64)
      @inbounds begin
            # Angular momentum
            vecadd!(kpl.sv.L, kpl.f.T, delta_t)
            # Angular velocity in global frame
            A_mul_B!(kpl.aux.ω, kpl.md.J_inv, kpl.sv.L)
            # Skew symmetric matrix
            skew4!(kpl.aux.ωs, kpl.aux.ω)
            # Quaternions time derivative
            q_vel!(kpl.dv.q̇, kpl.sv.q, kpl.aux.ωs)
            # Update & normalize Euler parameters
            q!(kpl.sv.q, kpl.dv.q̇, delta_t)
            normalize!(kpl.sv.q);
            # Update Rotation matrices
            rotmat!(kpl.aux.R, kpl.sv.q)
            inverse!(kpl.aux.R_inv, kpl.aux.R)
            # Update mass moment of inertia inverse in global frame
            Jb_inv2world!(kpl)
      end
      return nothing
end
"""
    body_dynamics!(kpl::Kappale, delta_t::Float64)
Calculates body's translation and rotation, and updates position in global frame. Updates all inplace.
"""
function body_dynamics!(kpl::Kappale, delta_t::Float64)
      body_translation!(kpl, delta_t)
      body_rotation!(kpl, delta_t)
      return nothing
end
