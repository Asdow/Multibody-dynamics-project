# -----------------------------
function Mmatrix(body::RigidBody{T}) where {T}
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
function Mmatrix!(M, body::RigidBody{T}) where {T}
    for i in 1:3
        @inbounds M[i,i] = body.md.m
    end
    for i in 4:6
        j = i-3;
        @inbounds M[i,i] = body.md.J[j,j]
    end
    return nothing
end
function Wmatrix(body::RigidBody{T}) where {T}
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
function Wmatrix!(W, body::RigidBody{T}) where {T}
    for i in 1:3
        W[i,i] = body.md.m_inv
    end
    for i in 4:6
        j = i-3;
        @inbounds W[i,i] = body.md.J_inv[j,j]
    end
    return nothing
end


function testi(Qe, bodies)
    stp = 1000000
    @time for i in 1:stp
        Qe!(Qe, bodies)
    end
    @time for i in 1:stp
        Qe2!(Qe, bodies)
    end
    return nothing
end


function test()
    Rsys = init_RSys([testicube(ones(3)..., 15.0) for i in 1:5]);
    nb = length(Rsys);
    bodies = Rsys.bodies;
    # set location and orientation to random
    randomize!(bodies)

    window = init_screen();
    origo = origin();
    bodyvis = [body_vis(bodies[i]) for i in 1:length(bodies)];

    gl._view(origo, window, camera=:perspective);
    gl._view(bodyvis, window, camera=:perspective);

    sleep(1)
    gv = sa.SVector{3,Float64}(1.0,2.0,3.0);
    # EOMs
    qddot = sa.MVector{6*nb,Float64}(zeros(6*nb));
    Qe = sa.MVector{6*nb,Float64}(zeros(6*nb));
    W = sysW(bodies);
    while Rsys.t < simTime
        clear_forces!(bodies);
        # vecadd!(body.f.F, g, body.md.m)
        if Rsys.t < 1.0
            for i in 1:nb
                vecadd!(bodies[i].f.T, gv, 1.0)
            end
        end
        translation!(bodies)
        rotation!(bodies)
        #-------------------
        # alternative version
        # sysWJ!(W, bodies)
        # Qe!(Qe, bodies)
        # A_mul_B!(qddot, W, Qe)
        # qddot2bodies!(bodies, qddot)
        # semiEuler!(bodies, nb)

        body_vis!(bodyvis, bodies)
        Rsys.t += Δt;
        sleep(0.001);
    end
    return nothing
end

function tesiti(bodies, nb)
    qddot = sa.MVector{6*nb,Float64}(zeros(6*nb));
    Qe = sa.MVector{6*nb,Float64}(zeros(6*nb));
    W = sysW(bodies);
    # Ws = sa.MMatrix{size(W,1), size(W,2), Float64, size(W,1)*size(W,2)}(W);
    stp = 10000
    @time for i in 1:stp
        translation!(bodies)
        rotation!(bodies)
    end
    @time for l in 1:stp
        sysWJ!(W, bodies)
        Qe!(Qe, bodies)
        A_mul_B!(qddot, W, Qe)
        qddot2bodies!(bodies, qddot)
        semiEuler!(bodies, nb)
    end
    return nothing
end

function semiEuler!(body::RigidBody{T}) where {T}
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
function semiEuler!(bodies, nb=length(bodies))
    for i in 1:nb
        @inbounds semiEuler!(bodies[i])
    end
    return nothing
end
semiEuler!(Rsys::RBodySystem) = semiEuler!(Rsys.bodies, Rsys.nb)
