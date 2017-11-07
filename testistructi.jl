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
import Base.eltype
eltype(body::RigidBody2{T}) where {T} = T
eltype(T2::Type{RigidBody2{T}}) where {T} = T

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
        @inbounds M[i,i] = body.md.m
    end
    for i in 4:6
        j = i-3;
        @inbounds M[i,i] = body.md.J[j,j]
    end
    return M
end
function Mmatrix!(M, body::RigidBody2{T}) where {T}
    for i in 1:3
        @inbounds M[i,i] = body.md.m
    end
    for i in 4:6
        j = i-3;
        @inbounds M[i,i] = body.md.J[j,j]
    end
    return nothing
end
function Wmatrix(body::RigidBody2{T}) where {T}
    W = sa.MMatrix{6,6,T,36}(zeros(T,6,6))
    for i in 1:3
        W[i,i] = body.md.m_inv
    end
    for i in 4:6
        j = i-3;
        @inbounds W[i,i] = body.md.J_inv[j,j]
    end
    return W
end
function Wmatrix!(W, body::RigidBody2{T}) where {T}
    for i in 1:3
        W[i,i] = body.md.m_inv
    end
    for i in 4:6
        j = i-3;
        @inbounds W[i,i] = body.md.J_inv[j,j]
    end
    return nothing
end

function sysM(bodies::Array{RigidBody2{T},1}) where {T<:Real}
    nb = length(bodies)
    S = 6*nb;
    M = zeros(T,S,S)
    sysMf!(M, i, bodies[i].md.m , bodies[i].md.J, nb)
    return M
end
function sysM!(M, bodies::Array{RigidBody2{T},1}) where {T}
    nb = length(bodies)
    for i in 1:nb
        sysMf!(M, i, bodies[i].md.m , bodies[i].md.J)
    end
    return nothing
end
function sysW(bodies::Array{RigidBody2{T},1}) where {T}
    nb = length(bodies)
    S = 6*nb;
    M = zeros(T,S,S)
    for i in 1:nb
        sysMf!(M, i, bodies[i].md.m_inv , bodies[i].md.J_inv)
    end
    return M
end
function sysW!(M, bodies::Array{RigidBody2{T},1}) where {T}
    nb = length(bodies)
    for i in 1:nb
        sysMf!(M, i, bodies[i].md.m_inv , bodies[i].md.J_inv)
    end
    return nothing
end

function sysMf!(M, i, m , J, nb)
    for i in 1:nb
        j = 1 + (6*(i-1))
        M[j,j] = m
        M[j+1,j+1] = m
        M[j+2,j+2] = m
        kk = 4 + (6*(i-1))
        for l in 0:2, ll in 0:2
            o = kk+l
            oo = kk+ll
            M[o,oo] = J[1+l, 1+ll]
        end
    end
    return nothing
end

function testi(Qe, bodies)
    stp = 1000
    @time for i in 1:stp
        Qe!(Qe, bodies)
    end
    # @time for i in 1:stp
    #     translation2!(bodies)
    # end
    return nothing
end
function foo(vektori)
    @time normalize!(vektori)
    return nothing
end
function translation!(body::RigidBody2{T}) where {T}
    body.ddv.a[:] = body.f.F.*body.md.m_inv
    body.dv.ẋ[:] += body.ddv.a .* Δt
    body.sv.x[:] += body.dv.ẋ .* Δt
    return nothing
end
function translation!(bodies::Array{RigidBody2{T},1}) where {T}
    for i in 1:length(bodies)
        @inbounds translation!(bodies[i])
    end
    return nothing
end

function rotation!(body::RigidBody2{T}) where {T}
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
function rotation!(bodies::Array{RigidBody2{T},1}) where {T}
    for i in 1:length(bodies)
        @inbounds rotation!(bodies[i])
    end
    return nothing
end

function test()
    bodies = [testicube(ones(3)..., 15.0) for i in 1:3];
    nb = length(bodies)
    # set location and orientation to random
    randomize!(bodies)

    window = init_screen();
    origo = origin();
    bodyvis = [body_vis(bodies[i]) for i in 1:length(bodies)];

    gl._view(origo, window, camera=:perspective);
    gl._view(bodyvis, window, camera=:perspective);

    sleep(1)
    t = 0.0;
    gv = sa.SVector{3,Float64}(1.0,2.0,3.0);
    # EOMs
    qddot = sa.MVector{6*nb,Float64}(zeros(6*nb));
    Qe = sa.MVector{6*nb,Float64}(zeros(6*nb));
    W = sysW(bodies);
    while t < simTime
        clear_forces!(body);
        # vecadd!(body.f.F, g, body.md.m)
        if t < 1.0
            for i in 1:nb
                vecadd!(bodies[i].f.T, gv, 1.0)
            end
        end
        # translation!(bodies)
        # rotation!(bodies)
        #-------------------
        # alternative version
        sysW!(W, bodies)
        Qe!(Qe, bodies)
        A_mul_B!(qddot, W, Qe)
        qddot2bodies!(bodies, qddot)
        for i in 1:nb
            bodies[i].dv.ẋ[:] += bodies[i].ddv.a .* Δt
            bodies[i].sv.x[:] += bodies[i].dv.ẋ .* Δt
            bodies[i].dv.ω[:] += bodies[i].ddv.α .* Δt
            skew4!(bodies[i].aux.ωs, bodies[i].dv.ω)
            # Quaternions time derivative
            q_vel!(bodies[i].dv.q̇, bodies[i].sv.q, bodies[i].aux.ωs)
            # Update & normalize Euler parameters
            bodies[i].sv.q[:] += bodies[i].dv.q̇ .* Δt
            normalize!(bodies[i].sv.q);
            # Update Rotation matrices
            rotmat!(bodies[i].aux.R, bodies[i].sv.q)
            inverse!(bodies[i].aux.R_inv, bodies[i].aux.R)
            # Update mass moment of inertias in global frame
            Jb2world!(bodies[i])
            Jb_inv2world!(bodies[i])
        end
        body_vis!(bodyvis, bodies)
        t += Δt;
        sleep(0.001);
    end
    return nothing
end

function randomize!(bodies)
    for i in 1:length(bodies)
        bodies[i].sv.x[:] = rand(3).*3+rand(3);
        bodies[i].sv.q[:] = rand(4);
        normalize!(bodies[i].sv.q)
        rotmat!(bodies[i])
        inverse!(bodies[i].aux.R_inv, bodies[i].aux.R)
        Jb2world!(bodies[i])
        Jb_inv2world!(bodies[i])
    end
    return nothing
end

function Qe!(Qe, bodies)
    nb = length(bodies)
    for i in 1:nb
        j = 1 + (6*(i-1))
        for o in 0:2
            Qe[j+o] = bodies[i].f.F[o+1];
        end
        # kk = (4:6) + (6*(i-1))
        # Qe[kk] = bodies[i].f.T - cross(bodies[i].dv.ω, (bodies[i].md.J*bodies[i].dv.ω))
        # Torque + Coriolis force
        kkk = 4 + (6*(i-1))
        A_mul_B!(bodies[i].ddv.a, bodies[i].md.J,bodies[i].dv.ω)
        Qe[kkk] = -( bodies[i].dv.ω[2]*bodies[i].ddv.a[3] - bodies[i].dv.ω[3]*bodies[i].ddv.a[2] )
        Qe[kkk+1] = -( bodies[i].dv.ω[3]*bodies[i].ddv.a[1] - bodies[i].dv.ω[1]*bodies[i].ddv.a[3] )
        Qe[kkk+2] = -( bodies[i].dv.ω[1]*bodies[i].ddv.a[2] - bodies[i].dv.ω[2]*bodies[i].ddv.a[1] )
        for o in 0:2
            Qe[kkk+o] += bodies[i].f.T[o+1]
        end
    end
    return nothing
end

function qddot2bodies!(bodies, qddot)
    for i in 1:length(bodies)
        j = 1 + (6*(i-1))
        kk = 4 + (6*(i-1))
        for o in 0:2
            bodies[i].ddv.a[o+1] = qddot[j+o]
            bodies[i].ddv.α[o+1] = qddot[kk+o]
        end
    end
    return nothing
end

function tesiti(bodies)
    qddot = sa.MVector{6*nb,Float64}(zeros(6*nb));
    Qe = sa.MVector{6*nb,Float64}(zeros(6*nb));
    W = sysW(bodies);
    stp = 1
    @time for i in 1:stp
        translation!(bodies)
        rotation!(bodies)
    end
    @time for l in 1:stp
        sysW!(W, bodies)
        Qe!(Qe, bodies)
        A_mul_B!(qddot, W, Qe)
        qddot2bodies!(bodies, qddot)
        for i in 1:nb
            @inbounds semiEuler(bodies[i])
        end
    end
    return nothing
end

function semiEuler(body)
    body.dv.ẋ[:] += body.ddv.a .* Δt
    body.sv.x[:] += body.dv.ẋ .* Δt
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
