# Datatypes used for the simulation
######################################
# Abstraktityyppi, jota käytetään tyyppien määrittelyssä
abstract type Kappale end
#################################################################################
# Datatyyppien määrittely Kappaleen konkreettisille tyypeille
mutable struct point3D{T<:Real} <: sa.FieldVector{3, T}
    x::T
    y::T
    z::T
end
mutable struct PMCoeff{T<:Real} <: sa.FieldVector{3, T}
      k::T # Contact stiffness
      c::T # Contact damping
      ki::T # Contact integration stiffness
end
struct AABB{T<:Real}
      min::point3D{T}
      max::point3D{T}
end
struct Shape{T<:Real}
      # Body shape described as a HomogenousMesh in local coordinates
      mesh::gt.HomogenousMesh{gt.Point{3,Float32},gt.Face{3,gt.OffsetInteger{-1,UInt32}},gt.Normal{3,Float32},Void,Void,Void,Void}
      # Global space AABB
      AABB::AABB{T}
      # Local space AABB
      AABBb::AABB{T}
end
function init_shape(x::T, y::T, z::T) where {T<:Real}
      body = cube(x, y, z)
      bb = AABB(point3D(zeros(T,3)), point3D(zeros(T,3)))
      sh = Shape(body, bb, deepcopy(bb) )
end
struct MassData{T<:Real}
      m::T # Mass
      m_inv::T
      Jb::sa.MArray{Tuple{3,3}, T, 2, 9} # Inertial matrix in local coordinates
      Jb_inv::sa.MArray{Tuple{3,3}, T, 2, 9}
      J::sa.MArray{Tuple{3,3}, T, 2, 9} # Inertial matrix in global coordinates
      J_inv::sa.MArray{Tuple{3,3}, T, 2, 9}
end
function init_MassData(m::T) where {T<:Real}
      Jb = sa.MArray{Tuple{3,3},T}(zeros(T,3,3))
      # Handle 0 mass special case. Creates a body with inverted mass == 0. Won't move anywhere.
      if iszero(m) == true
            MD = MassData(m, m, Jb, deepcopy(Jb), deepcopy(Jb), deepcopy(Jb))
      else
            MD = MassData(m, inv(m), Jb, deepcopy(Jb), deepcopy(Jb), deepcopy(Jb))
      end
      return MD
end
struct StateVariables{T<:Real}
      x::point3D{T} # Reference point location in global coordinates
      q::mQuaternion{T} # Euler parameters
end
function init_StateVariables(T)
      q = mQuaternion(zeros(T,4))
      q.s = 1.0;
      sv = StateVariables( point3D(zeros(T,3)), q )
end
struct Derivatives{T<:Real}
      ẋ::sa.MVector{3,T} # Reference point velocity in global coordinates
      q̇::mQuaternion{T} # Euler parameters time derivative
      ω::sa.MVector{3,T} # Body angular velocity in global coordinates
end
function init_Derivatives(T)
    vec = sa.MVector{3, T}(zeros(T,3))
    dv = Derivatives( vec, mQuaternion(zeros(T,4)), deepcopy(vec) )
end
struct SecDerivatives{T<:Real}
    a::sa.MVector{3,T} # Reference point acc in global coordinates
    α::sa.MVector{3,T} # ang. acc.
end
function init_SecDerivatives(T)
    vec = sa.MVector{3, T}(zeros(T,3))
    dv = SecDerivatives( vec, deepcopy(vec) )
end
struct AuxValues{T<:Real}
      R::sa.MArray{Tuple{3,3}, T, 2, 9} # Rotation matrix. Local -> global
      R_inv::sa.MArray{Tuple{3,3}, T, 2, 9} # Rotation matrix. Global -> local
      ωs::sa.MArray{Tuple{4,4}, T, 2, 16} # Skew symmetric matrix
end
function init_Aux(T)
      Jb = sa.MArray{Tuple{3,3},T}(zeros(T,3,3))
      Jb2 = sa.MArray{Tuple{4,4},T}(zeros(T,4,4))
      aux = AuxValues( Jb, deepcopy(Jb), Jb2 )
end
struct Forces{T<:Real}
      F::sa.MVector{3,T} # Forces effecting the body in global coordinates
      T::sa.MVector{3,T} # Moments effecting the body in global coordinates
end
function init_Forces(T)
      vec = sa.MVector{3, T}(zeros(T,3))
      forces = Forces( vec, deepcopy(vec) )
end
struct PenaltyMethod{T<:Real}
      coeff::PMCoeff{T}
      max_dyi::sa.MVector{100,T} # Max penetration error
end
function init_PenMethod(T)
      pm = PenaltyMethod( PMCoeff(zeros(T,3)), sa.MVector{100, T}(zeros(T,100)) )
end
struct LuGre{T<:Real}
      Fs::sa.MVector{3,T} # Static friction force
      Fk::sa.MVector{3,T} # Kinetic friction force
      Ts::sa.MVector{3,T} # Static friction moment
      Tk::sa.MVector{3,T} # Kinetic friction moment
      μs::T # Static friction coefficient
      μk::T # Kinetic friction coefficient
end
function init_Lugre(mus::T, muk::T) where {T<:Real}
      vec = sa.MVector{3, T}(zeros(T,3))
      lugre = LuGre( vec, deepcopy(vec), deepcopy(vec), deepcopy(vec), mus, muk )
end
#################################################################################
# Kappaletyyppien määrittely
struct RigidBody{T<:Real} <: Kappale
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

mutable struct RBodySystem{T<:Real}
    bodies::Array{RigidBody{T},1} # A list of bodies
    nb::Int64 # Number of bodies
    t::Float64 # time
end
function init_RSys(bodies::Array{RigidBody{T},1}) where {T}
    Rsys = RBodySystem(bodies, length(bodies), 0.0)
end
import Base.length
length(a::RBodySystem) = a.nb
