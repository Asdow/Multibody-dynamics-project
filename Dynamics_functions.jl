# Functions for the rigid body dynamics
"""
    clear_forces!(kpl::T) where {T<:Kappale}
Clears body's forces.
"""
function clear_forces!(kpl::T) where {T<:Kappale}
    fill!(kpl.f.F, 0.0);
    fill!(kpl.f.T, 0.0);
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
    eval_forces!(kpl::T) where {T<:Kappale}
Evaluates RBodySystem's forces and adds them to bodies' force accumulators.
"""
function eval_forces!(kpl::T) where {T<:Kappale}
    clear_forces!(kpl)
    gravity!(kpl)
    return nothing
end
function bodies::Array{T,1}, nb=length(bodies)) where {T<:Kappale}
    for i in 1:nb
        eval_forces!(bodies[i])
    end
    return nothing
end
eval_forces!(Rsys::RBodySystem{T}) where {T} = eval_forces!(Rsys.bodies, Rsys.nb)


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
    translation!(body::RigidBody{T}) where {T}
Calculates body's linear acceleration, velocity and reference point global location. Updates all three inplace.
"""
function translation!(body::RigidBody{T}) where {T}
    body.ddv.a[:] = body.f.F.*body.md.m_inv
    body.dv.ẋ[:] += body.ddv.a .* Δt
    body.sv.x[:] += body.dv.ẋ .* Δt
    return nothing
end
function translation!(bodies::Array{RigidBody{T},1}, nb=length(bodies)) where {T}
    for i in 1:nb
        @inbounds translation!(bodies[i])
    end
    return nothing
end
translation!(Rsys::RBodySystem) = translation!(Rsys.bodies, Rsys.nb)
"""
    rotation!(body::RigidBody{T}) where {T}
Calculates body's angular acceleration, velocity, euler parameters rate of change, and reference point global orientation. Updates all inplace.
"""
function rotation!(body::RigidBody{T}) where {T}
    body.ddv.α[:] = body.md.J_inv*(body.f.T - cross(body.dv.ω, (body.md.J*body.dv.ω)))
    body.dv.ω[:] += body.ddv.α .* Δt
    skew4!(body.aux.ωs, body.dv.ω)
    # Quaternions time derivative
    q_vel!(body.dv.q̇, body.sv.q, body.aux.ωs)
    # Update & normalize Euler parameters
    body.sv.q[:] += body.dv.q̇ .* Δt
    normalize!(body.sv.q);
    # Update Rotation matrices
    rotmat!(body.aux.R, body.sv.q)
    inverse!(body.aux.R_inv, body.aux.R)
    # Update mass moment of inertias in global frame
    Jb2world!(body)
    Jb_inv2world!(body)
    return nothing
end
function rotation!(bodies::Array{RigidBody{T},1}, nb=length(bodies)) where {T}
    for i in 1:nb
        @inbounds rotation!(bodies[i])
    end
    return nothing
end
rotation!(Rsys::RBodySystem) = rotation!(Rsys.bodies, Rsys.nb)
"""
    dynamics!(kpl::Kappale)
Calculates body's translation and rotation, and updates position in global frame. Updates all inplace.
"""
function dynamics!(kpl::Kappale)
      translation!(kpl)
      rotation!(kpl)
      return nothing
end
function dynamics!(bodies::Array{RigidBody{T},1}, nb=length(bodies)) where {T}
    for i in 1:nb
        @inbounds dynamics!(bodies[i])
    end
    return nothing
end
dynamics!(Rsys::RBodySystem) = dynamics!(Rsys.bodies, Rsys.nb)
"""
    gravity!(body::Kappale)
Calculates a body's gravity force and adds it to its force accumulator.
"""
function gravity!(body::Kappale)
    body.f.F[:] += g.*body.md.m
    return nothing
end
function gravity!(bodies::Array{RigidBody{T},1}, nb=length(bodies)) where {T}
    for i in 1:nb
        @inbounds gravity!(bodies[i])
    end
    return nothing
end
gravity!(Rsys::RBodySystem) = gravity!(Rsys.bodies, Rsys.nb)
