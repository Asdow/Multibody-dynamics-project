struct SV{T<:Real}
      x::point3D{T} # Reference point location in global coordinates
      q::mQuaternion{T} # Euler parameters
end
function init_SV(T)
      q = mQuaternion(zeros(T,4))
      q.s = 1.0;
      sv = SV( point3D(zeros(T,3)), q)
end
struct DV{T<:Real}
      ẋ::sa.MVector{3,T} # Reference point velocity in global coordinates
      q̇::mQuaternion{T} # Euler parameters time derivative
      ω::sa.MVector{3,T} # ang. vel.
end
function init_DV(T)
    vec = sa.MVector{3, T}(zeros(T,3))
    dv = DV( vec, mQuaternion(zeros(T,4)), deepcopy(vec) )
end
struct DDV{T<:Real}
      a::sa.MVector{3,T} # Reference point acc in global coordinates
      α::sa.MVector{3,T} # ang. acc.
end
function init_DDV(T)
    vec = sa.MVector{3, T}(zeros(T,3))
    dv = DDV( vec, deepcopy(vec) )
end
struct RigidBody2{T<:Real} <: Kappale
      sh::Shape{T}
      md::MassData{T}
      sv::SV{T}
      dv::DV{T}
      ddv::DDV{T}
      aux::AuxValues{T}
      f::Forces{T}
end

function testicube(x::T, y::T, z::T, m::T) where {T<:Real}
      sh = init_shape(x, y, z)
      md = init_MassData(m)
      sv = init_SV(T)
      dv = init_DV(T)
      ddv = init_DDV(T)
      aux = init_Aux(T)
      f = init_Voimat(T)
      body = RigidBody2( sh, md, sv, dv, ddv, aux, f )
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

# -----------------------------
function Mmatrix(body::RigidBody2{T}) where {T}
    M = sa.MMatrix{6,6,T,36}(zeros(T,6,6))
    for i in 1:3
        M[i,i] = body.md.m
    end
    M[4:6,4:6] = body.md.J
    return M
end
function Wmatrix(body::RigidBody2{T}) where {T}
    W = sa.MMatrix{6,6,T,36}(zeros(T,6,6))
    for i in 1:3
        W[i,i] = body.md.m_inv
    end
    W[4:6,4:6] = body.md.J_inv
    return W
end

function translation(body::RigidBody2{T}) where {T}
    body.ddv.a[:] = body.f.F.*body.md.m_inv
    body.dv.ẋ[:] += body.ddv.a .* delta_t
    body.sv.x[:] += body.dv.ẋ .* delta_t
    return nothing
end

function rotation(body::RigidBody2{T}) where {T}
    body.ddv.α[:] = body.md.J_inv*(body.f.T - cross(body.dv.ω, (body.md.J*body.dv.ω)))
    body.dv.ω[:] += body.ddv.α .* delta_t
    skew4!(body.aux.ωs, body.dv.ω)
    # Quaternions time derivative
    q_vel!(body.dv.q̇, body.sv.q, body.aux.ωs)
    # Update & normalize Euler parameters
    body.sv.q[:] += body.dv.q̇ .* delta_t
    normalize!(body.sv.q);
    # Update Rotation matrices
    rotmat!(body.aux.R, body.sv.q)
    inverse!(body.aux.R_inv, body.aux.R)
    # Update mass moment of inertias in global frame
    Jb2world!(body)
    Jb_inv2world!(body)
    return nothing
end

function test()
    body = testicube(ones(3)..., 15.0);
    # set location and orientation to random
    body.sv.x[:] = rand(3);
    body.sv.q[:] = rand(4);
    normalize!(body.sv.q)
    rotmat!(body)
    inverse!(body.aux.R_inv, body.aux.R)
    Jb2world!(body)
    Jb_inv2world!(body)

    window = init_screen();
    origo = origin();
    bodyvis = body_vis(body);

    gl._view(origo, window, camera=:perspective);
    gl._view(bodyvis, window, camera=:perspective);

    sleep(1)
    t = 0.0;
    gv = sa.SVector{3,Float64}(1.0,2.0,3.0);
    qddot = sa.MVector{6,Float64}(zeros(6));
    Qe = sa.MVector{6,Float64}(zeros(6));
    while t < simTime
        clear_forces!(body);
        # vecadd!(body.f.F, g, body.md.m)
        # if t < 1.0
        #     vecadd!(body.f.T, gv, 1.0)
        # end
        # translation(body)
        # rotation(body)
        M = Mmatrix(body)
        Qe[1:3] = body.f.F;
        Qe[4:6] = body.f.T - cross(body.dv.ω, (body.md.J*body.dv.ω))
        qddot[:] = inv(M)*Qe
        body.ddv.a[:] = qddot[1:3]
        body.ddv.α[:] = qddot[4:6]

        body.dv.ẋ[:] += body.ddv.a .* delta_t
        body.sv.x[:] += body.dv.ẋ .* delta_t
        body.dv.ω[:] += body.ddv.α .* delta_t
        skew4!(body.aux.ωs, body.dv.ω)
        # Quaternions time derivative
        q_vel!(body.dv.q̇, body.sv.q, body.aux.ωs)
        # Update & normalize Euler parameters
        body.sv.q[:] += body.dv.q̇ .* delta_t
        normalize!(body.sv.q);
        # Update Rotation matrices
        rotmat!(body.aux.R, body.sv.q)
        inverse!(body.aux.R_inv, body.aux.R)
        # Update mass moment of inertias in global frame
        Jb2world!(body)
        Jb_inv2world!(body)

        gl.set_arg!(bodyvis, :model, transformation(body))
        t += delta_t;
        sleep(0.001);
    end
    return nothing
end
