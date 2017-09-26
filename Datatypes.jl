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
      body::gt.HomogenousMesh{gt.Point{3,Float32},gt.Face{3,gt.OffsetInteger{-1,UInt32}},gt.Normal{3,Float32},Void,Void,Void,Void}
      # Global space AABB
      AABB::AABB{T}
      # Local space AABB
      AABBb::AABB{T}
end
struct MassData{T<:Real}
      m::T # Mass
      m_inv::T
      Jb::sa.MArray{Tuple{3,3}, T, 2, 9} # Inertial matrix in local coordinates
      Jb_inv::sa.MArray{Tuple{3,3}, T, 2, 9}
      J::sa.MArray{Tuple{3,3}, T, 2, 9} # Inertial matrix in global coordinates
      J_inv::sa.MArray{Tuple{3,3}, T, 2, 9}
end
struct StateVariables{T<:Real}
      x::point3D{T} # Reference point location in global coordinates
      q::mQuaternion{T} # Euler parameters
      P::sa.MVector{3,T} # Linear momentum in global coordinates
      L::sa.MVector{3,T} # Angular momentum in global coordinates
end
struct Derivates{T<:Real}
      ẋ::point3D{T} # Reference point velocity in global coordinates
      q̇::mQuaternion{T} # Euler parameters time derivative
end
struct AuxValues{T<:Real}
      R::sa.MArray{Tuple{3,3}, T, 2, 9} # Rotation matrix. Local -> global
      R_inv::sa.MArray{Tuple{3,3}, T, 2, 9} # Rotation matrix. Global -> local
      ωb::sa.MVector{3,T} # Body angular velocity in local coordinates
      ω::sa.MVector{3,T} # Body angular velocity in global coordinates
      ωs::sa.MArray{Tuple{4,4}, T, 2, 16} # Skew symmetric matrix
end
struct Forces{T<:Real}
      Fb::sa.MVector{3,T} # Forces effecting the body in local coordinates
      F::sa.MVector{3,T} # Forces effecting the body in global coordinates
      Tb::sa.MVector{3,T} # Moments effecting the body in local coordinates
      T::sa.MVector{3,T} # Moments effecting the body in global coordinates
end
struct PenaltyMethod{T<:Real}
      coeff::PMCoeff{T}
      max_dyi::sa.MVector{100,T} # Max penetration error
end
struct LuGre{T<:Real}
      Fs::sa.MVector{3,T} # Static friction force
      Fk::sa.MVector{3,T} # Kinetic friction force
      Ts::sa.MVector{3,T} # Static friction moment
      Tk::sa.MVector{3,T} # Kinetic friction moment
      μs::T # Static friction coefficient
      μk::T # Kinetic friction coefficient
end
#################################################################################
# Kappaletyyppien määrittely
struct RigidBody <: Kappale
      sh::Shape
      md::MassData
      sv::StateVariables
      dv::Derivates
      aux::AuxValues
      f::Forces
      knt::PenaltyMethod
      kd::LuGre
end
######################################
# Init function for Datatypes
######################################
# Shapedata
function init_shape(x::T, y::T, z::T) where {T<:Real}
      body = cube(x, y, z)
      bb = AABB(point3D(zeros(3)), point3D(zeros(3)))
      sh = Shape(body, bb, deepcopy(bb) )
end
# MassData
function init_MassData(m::Float64)
      Jb = sa.MArray{Tuple{3,3},Float64}(zeros(3,3))
      # Handle 0 mass special case. Creates a body with inverted mass == 0. Won't move anywhere.
      if iszero(m) == true
            MD = MassData(m, m, Jb, deepcopy(Jb), deepcopy(Jb), deepcopy(Jb))
      else
            MD = MassData(m, inv(m), Jb, deepcopy(Jb), deepcopy(Jb), deepcopy(Jb))
      end
      return MD
end
# StateVariables
function init_StateVariables()
      vec = sa.MVector{3, Float64}(zeros(3))
      q = mQuaternion(zeros(4))
      q.s = 1.0;
      sv = StateVariables( point3D(zeros(3)), q, vec, deepcopy(vec) )
end
# Derivates
function init_Derivates()
      dv = Derivates( point3D(zeros(3)), mQuaternion(zeros(4)) )
end
# AuxValues
function init_Aux()
      Jb = sa.MArray{Tuple{3,3},Float64}(zeros(3,3))
      Jb2 = sa.MArray{Tuple{4,4},Float64}(zeros(4,4))
      vec = sa.MVector{3, Float64}(zeros(3))
      aux = AuxValues( Jb, deepcopy(Jb), vec, deepcopy(vec), Jb2 )
end
# Forces
function init_Voimat()
      vec = sa.MVector{3, Float64}(zeros(3))
      forces = Forces( vec, deepcopy(vec), deepcopy(vec), deepcopy(vec) )
end
# PenaltyMethod
function init_PenMethod()
      pm = PenaltyMethod( PMCoeff(zeros(3)), sa.MVector{100, Float64}(zeros(100)) )
end
# LuGre
function init_Lugre(mus, muk)
      vec = sa.MVector{3, Float64}(zeros(3))
      lugre = LuGre( vec, deepcopy(vec), deepcopy(vec), deepcopy(vec), mus, muk )
end
