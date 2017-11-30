struct RigidBody{T} <: Kappale{T}
      sh::Shape{T}
      md::MassData{T}
      sv::StateVariables{T}
      dv::Derivatives{T}
      ddv::SecDerivatives{T}
      aux::AuxValues{T}
      f::Forces{T}
      knt::PenaltyMethod{T}
      kd::LuGre{T}
end
import Base.eltype
eltype(body::RigidBody{T}) where {T} = T
eltype(T2::Type{RigidBody{T}}) where {T} = T

getWorldLocation(body::Kappale{T}) where {T} = getx(body.sv)
getWorldOrientation(body::Kappale{T}) where {T} = getq(body.sv)

function randomize!(body::T) where {T<:Kappale}
    x = getWorldLocation(body)
    q = getWorldOrientation(body)

    x[:] = rand(3).*3+rand(3);
    q[:] = rand(4);
    normalize!(q)
    rotmat!(body)
    inverse!(getRotmatinv(body), getRotmat(body))
    Inertias2world!(body)
    return nothing
end
function randomize!(bodies, nb=length(bodies))
    for i in 1:nb
        randomize!(bodies[i])
    end
    return nothing
end

getShape(body::RigidBody{T}) where {T} = body.sh
getCollShapes(body::RigidBody{T}) where {T} = getShape(body).coll

function setCollShapes!(body::RigidBody{T}, cs::Collshape{T}) where {T}
    sh = getShape(body)
    setCollShapes!(sh, cs)
    return nothing
end
function setCollShapes!(sh::Shape{T}, cs::Collshape{T}) where {T}
    sh.coll = cs
    return nothing
end
