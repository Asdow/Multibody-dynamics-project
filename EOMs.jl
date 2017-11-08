struct EOM{T<:Real}
      W::Array{T,2} # Inverse of system Mass Matrix
      Qe::Array{T,1} # Vector of external forces
      ü::Array{T,1} # Vector of gen. accelerations
end
function init_EOM(Rsys::RBodySystem{T}) where {T}
    # S = nb*6
    # L = S*S
    # W = sa.MMatrix{S, S, T, L}(zeros(T, S,S));
    # Qe = sa.MVector{S,T}(zeros(T,S));
    # ü = deepcopy(Qe);
    W = sysW(Rsys)
    Qe = sysQe(Rsys)
    ü = zeros(Qe)
    EOMs = EOM(W, Qe, ü)
end
#-----------------------------------------------------------
function sysM(bodies::Array{RigidBody{T},1}, nb=length(bodies)) where {T<:Real}
    # nb = length(bodies)
    S = 6*nb;
    M = zeros(T,S,S)
    for i in 1:nb
        sysMf!(M, bodies[i].md.m , bodies[i].md.J, i)
    end
    return M
end
sysM(Rsys::RBodySystem) = sysM(Rsys.bodies, Rsys.nb)
function sysM!(M, bodies::Array{RigidBody{T},1}, nb=length(bodies)) where {T}
    # nb = length(bodies)
    for i in 1:nb
        sysMf!(M, bodies[i].md.m, bodies[i].md.J, i)
    end
    return nothing
end
sysM!(M,Rsys::RBodySystem) = sysM!(M,Rsys.bodies, Rsys.nb)
function sysW(bodies::Array{RigidBody{T},1}, nb=length(bodies)) where {T}
    # nb = length(bodies)
    S = 6*nb;
    M = zeros(T,S,S)
    for i in 1:nb
        sysMf!(M, bodies[i].md.m_inv, bodies[i].md.J_inv, i)
    end
    # Ms = sa.MMatrix{size(M,1), size(M,2), Float64, size(M,1)*size(M,2)}(M);
    return M
end
sysW(Rsys::RBodySystem) = sysW(Rsys.bodies, Rsys.nb)
function sysW!(M, bodies::Array{RigidBody{T},1}, nb=length(bodies)) where {T}
    # nb = length(bodies)
    for i in 1:nb
        sysMf!(M, bodies[i].md.m_inv, bodies[i].md.J_inv, i)
    end
    return nothing
end
sysW!(W,Rsys::RBodySystem) = sysW!(W, Rsys.bodies, Rsys.nb)
sysW!(EOMs::EOM, Rsys::RBodySystem) = sysW!(EOMs.W,Rsys)

function sysWJ!(M, bodies::Array{RigidBody{T},1}, nb=length(bodies)) where {T}
    # nb = length(bodies)
    for i in 1:nb
        sysMfJ!(M, bodies[i].md.J_inv, i)
    end
    return nothing
end
sysWJ!(W,Rsys::RBodySystem) = sysWJ!(W, Rsys.bodies, Rsys.nb)
sysWJ!(EOMs::EOM, Rsys::RBodySystem) = sysWJ!(EOMs.W,Rsys)

function sysMf!(M, m, J, i)
    sysMfm!(M, m, i)
    sysMfJ!(M, J, i)
    return nothing
end
function sysMfm!(M, m, i)
    j = 1 + (6*(i-1))
    M[j,j] = m
    M[j+1,j+1] = m
    M[j+2,j+2] = m
    return nothing
end
function sysMfJ!(M, J, i)
    kk = 4 + (6*(i-1))
    for l in 0:2, ll in 0:2
        o = kk+l
        oo = kk+ll
        M[o,oo] = J[1+l, 1+ll]
    end
    return nothing
end
function sysQe(bodies::Array{RigidBody{T},1}, nb=length(bodies)) where {T}
    Qe = zeros(T, 6*nb)
    for i in 1:nb
        QeF!(Qe, bodies[i], i)
        QeT!(Qe, bodies[i], i)
    end
    return Qe
end
sysQe(Rsys::RBodySystem) = sysQe(Rsys.bodies, Rsys.nb)
function sysQe!(Qe, bodies, nb=length(bodies))
    # nb = length(bodies)
    for i in 1:nb
        QeF!(Qe, bodies[i], i)
        QeT!(Qe, bodies[i], i)
    end
    return nothing
end
sysQe!(Qe, Rsys::RBodySystem) = sysQe!(Qe, Rsys.bodies, Rsys.nb)

function QeF!(Qe, body, i)
    j = 1 + (6*(i-1))
    for o in 0:2
        Qe[j+o] = body.f.F[o+1];
    end
    return nothing
end
function QeT!(Qe, body, i)
    # kk = (4:6) + (6*(i-1))
    # Qe[kk] = bodies[i].f.T - cross(bodies[i].dv.ω, (bodies[i].md.J*bodies[i].dv.ω))
    # Torque + Coriolis force
    kkk = 4 + (6*(i-1))
    A_mul_B!(body.ddv.a, body.md.J,body.dv.ω)
    Qe[kkk] = -( body.dv.ω[2]*body.ddv.a[3] - body.dv.ω[3]*body.ddv.a[2] )
    Qe[kkk+1] = -( body.dv.ω[3]*body.ddv.a[1] - body.dv.ω[1]*body.ddv.a[3] )
    Qe[kkk+2] = -( body.dv.ω[1]*body.ddv.a[2] - body.dv.ω[2]*body.ddv.a[1] )
    for o in 0:2
        Qe[kkk+o] += body.f.T[o+1]
    end
    return nothing
end

function qddot2bodies!(bodies, qddot, nb=length(bodies))
    for i in 1:nb
        j = 1 + (6*(i-1))
        kk = 4 + (6*(i-1))
        for o in 0:2
            bodies[i].ddv.a[o+1] = qddot[j+o]
            bodies[i].ddv.α[o+1] = qddot[kk+o]
        end
    end
    return nothing
end
