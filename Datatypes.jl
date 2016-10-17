# Datatypes used for the simulation
######################################
using StaticArrays
using ArrayFire
# Abstraktityyppi, jota käytetään tyyppien määrittelyssä
abstract Kappale
#################################################################################
# Datatyyppien määrittely Kappaleen konkreettisille tyypeille
type PointData
      n::Integer              # Kappaleen pisteiden lukumäärä
      k::Integer              # Kappaleen pisteiden lukumäärä
      w::AbstractFloat            # Kappaleen pituus
end

type Worldcoords
      # Kappaleen asema ja orientiaatio globaalissa koordinaatistossa
      XYZ::AbstractArray  # Kappaleen pisteiden asema globaalissa koordinaatistossa
      X::AbstractArray  # Kappaleen pisteiden X-asema globaalissa koordinaatistossa
      Y::AbstractArray  # Kappaleen pisteiden Y-asema globaalissa koordinaatistossa
      Z::AbstractArray  # Kappaleen pisteiden Z-asema globaalissa koordinaatistossa
end

type Bodycoords
      # Kappaleen muoto kuvataan lokaalissa koordinaatistossa
      R::AbstractArray      # Kappaleen muoto lokaalissa sylinterikoordinaatistossa, säde
      Θ::AbstractArray   # Kappaleen muoto lokaalissa sylinterikoordinaatistossa, kulma
      XYZ::AbstractArray  # Kappaleen pisteiden asema lokaalissa koordinaatistossa
end

type MassData
      ρ::AbstractFloat      # Kappaleen tiheys
      volume::AbstractFloat   # Kappaleen tilavuus
      massa::AbstractFloat    # Kappaleen massa
end

type Inertias
      body::AbstractArray # Kappaleen hitausmomentti lokaalissa koordinaatistossa
      body_inv::AbstractArray
      world::AbstractArray # Kappaleen hitausmomentti globaalissa koordinaatistossa
      world_inv::AbstractArray
end

type StateVariables
      x::AbstractArray #StaticArrays.MVector{3,Float64}        # Massakeskipisteen asema globaalissa koordinaatistossa
      q::AbstractArray #StaticArrays.MVector{4,Float64}        # Eulerin parametrit
      P::AbstractArray #StaticArrays.MVector{3,Float64}       # Lineaarinen liikemäärä globaalissa koordinaatistossa
      L::AbstractArray #StaticArrays.MVector{3,Float64}        # Angulaarinen liikemäärä
end

type Derivaatat
      ẋ::AbstractArray              # Massakeskipisteen lineaarinen nopeus globaalissa koordinaatistossa
      ẍ::AbstractArray              # Massakeskipisteen lineaarinen kiihtyvyys globaalissa koordinaatistossa
      q̇::AbstractArray             # Eulerin parametrien muutosnopeus
      q̈::AbstractArray             # Eulerin parametrien kiihtyvyys
end

type ApuArvot
      Rot_mat::AbstractArray         # Rotaatiomatriisi
      Rot_mat_inv::AbstractArray     # Rotaatiomatriisin käänteismatriisi
      ω_body::AbstractArray     # Kappaleen kulmanopeus lokaalissa koordinaatistossa
      ω_world::AbstractArray    # Kappaleen kulmanopeus globaalissa koordinaatistossa
      ω_tilde::AbstractArray     # 3x3 skew symmetric omega_tilde matriisi
      ω_tilde4::AbstractArray    # 4x4 skew symmetric omega_tilde matriisi
      E::AbstractArray               # 4x4 E matriisi
      E_1::AbstractArray             # 3x4 E matriisi
end

type Voimat
      force_body::AbstractArray         # Kappaleeseen kohdistuvat voimat lokaalissa koordinaatistossa
      force_world::AbstractArray         # Kappaleeseen kohdistuvat voimat globaalissa koordinaatistossa
      torque_body::AbstractArray   # Kappaleeseen kohdistuvat väännöt lokaalissa koordinaatistossa
      torque_world::AbstractArray  # Kappaleeseen kohdistuvat väännöt globaalissa koordinaatistossa
end

type KontaktiData
      k::AbstractFloat      # Kappaleen jousivakio kontaktivoimien laskemiseen
      c::AbstractFloat      # Kappaleen vaimennusvakio kontaktivoimien laskemiseen
      vv::AbstractArray # Matriisi[{Vector3}, n+1,k], joka säilyttää kappaleen kontaktista aiheutuvat voimavektorit pöllin pisteille
      vv_perpendicular::AbstractArray # Matriisi, joka säilyttää voimavektoreita kohtisuoraan olevat voimavektorit
      torque_vektorit::AbstractArray # Matriisi, joka säilyttää voimavektoreista aiheutuvat pölliin kohdistuvat momentit.
      F_penalty::AbstractArray # Kontakteista johtuvat rangaistusvoimat
      painuma_max::AbstractArray # Kontaktien painuma
      painuma_max_edellinen::AbstractArray # Kontaktien painuma edellisellä aika-askeleella
end

type BoundingBox
      # Bounding box matriisit kappaleelle.
      x::AbstractArray # Vektori, jossa 1. solu on min-arvo ja 2. max-arvo
      y::AbstractArray
      z::AbstractArray
      box::AbstractArray # Matriisi[2,3], 1. sarake on x:n min ja max-arvot, 2. y:n  ja 3. z:n.
end

type KitkaData
      F_friction_static::AbstractFloat   # Suurin staattinen kitkavoima
      F_friction_kinetic::AbstractFloat  # Suurin kineettinen kitkavoima
      T_friction_static::AbstractFloat   # Suurin staattinen kitkamomentti
      T_friction_kinetic::AbstractFloat  # Suurin kineettinen kitkamomentti
      μ_s::AbstractFloat            # Staattinen kitkakerroin
      μ_k::AbstractFloat            # Kineettinen kitkakerroin
      Slipping::Integer         # Int-arvo joka kertoo luistaako rulla vai ei
      F_torque::AbstractFloat      # Tukirullan väännöstä johtuva voima, joka välittyy pölliin kitkan avulla
      T_friction::AbstractFloat     # Kitkavoimasta johtuva, pölliin välitetty vääntömomentti
end

#################################################################################
# Kappaletyyppien määrittely
type RigidBody <: Kappale
      # Pistedata
      pd::PointData
      #####
      # Kappaleen koordinaattidata
      body::Bodycoords
      world::Worldcoords
      #####
      # Kappaleen massa ja tilavuusdata
      md::MassData
      #####
      # Hitausmomentit
      J::Inertias
      #####
      # State variables
      sv::StateVariables
      #####
      # Derivaatat
      dv::Derivaatat
      #####
      # Apuarvot
      aux::ApuArvot
      #####
      # Lasketut arvot
      f::Voimat
      #####
      # KontaktiData
      knt::KontaktiData
      #####
      # Kitkadata
      kd::KitkaData
      #####
      # bounbing box
      bb::BoundingBox
      #####
end

#################################################################################
# Init funktiot Datatyypeille
#####
# PointData
function func_init_pd(n::Integer, k::Integer, w::AbstractFloat)
      Pointdatatemp = PointData(n, k, w)
end
######################################
# Worldcoords
function func_init_Worldcoords(n::Integer, k::Integer, T::Type)
      tempvector3 = zeros(T,3)
      temparray2 = zeros(T, n+1,k)
      temparray3 = Array(Array{T},n+1,k)
      @simd for j in range(1,k)
          @simd for i in range(1,n+1)
              @inbounds temparray3[i,j] = deepcopy(tempvector3)
          end
      end

      Worldcoordstemp = Worldcoords(temparray3, deepcopy(temparray2), deepcopy(temparray2), deepcopy(temparray2) )
end
function func_init_Worldcoords(n::Integer, k::Integer, T::Type{StaticArray})
      tempvector3 = @MVector(zeros(3))
      temparray2 = @MMatrix(zeros(n+1,k))
      temparray3 = Array(StaticArrays.MVector{3,Float64},n+1,k)
      @simd for j in range(1,k)
          @simd for i in range(1,n+1)
              @inbounds temparray3[i,j] = deepcopy(tempvector3)
          end
      end

      Worldcoordstemp = Worldcoords(temparray3, deepcopy(temparray2), deepcopy(temparray2), deepcopy(temparray2) )
end
function func_init_Worldcoords(n::Integer, k::Integer, T::Type{AFArray})
      temparray2 = zeros(n+1,k)
      temparray3 = Array(ArrayFire.AFArray{Float64},n+1,k)
      @simd for j in range(1,k)
          @simd for i in range(1,n+1)
              @inbounds temparray3[i,j] = AFArray(zeros(3))
          end
      end

      Worldcoordstemp = Worldcoords(temparray3, AFArray(temparray2), AFArray(temparray2), AFArray(temparray2) )
end
######################################
# Bodycoords
function func_init_Bodycoords(n::Integer, k::Integer, T::Type)
      tempvector3 = zeros(T,3)
      tempvector = zeros(T, n+1)
      temparray2 = zeros(T, n+1,k)
      temparray3 = Array(Array{T},n+1,k)
      @simd for j in range(1,k)
          @simd for i in range(1,n+1)
              @inbounds temparray3[i,j] = deepcopy(tempvector3)
          end
      end

      Bodycoordstemp = Bodycoords( temparray2, tempvector, temparray3 )
end
function func_init_Bodycoords(n::Integer, k::Integer, T::Type{StaticArray})
      tempvector3 = @MVector(zeros(3))
      temparray2 = @MMatrix(zeros(n+1,k))
      temparray3 = Array(StaticArrays.MVector{3,Float64},n+1,k)
      tempvector = @MVector(zeros(n+1))
      @simd for j in range(1,k)
          @simd for i in range(1,n+1)
              @inbounds temparray3[i,j] = deepcopy(tempvector3)
          end
      end

      Bodycoordstemp = Bodycoords( temparray2, tempvector, temparray3 )
end
function func_init_Bodycoords(n::Integer, k::Integer, T::Type{AFArray})
      temparray2 = zeros(n+1,k)
      temparray3 = Array(ArrayFire.AFArray{Float64},n+1,k)
      tempvector = zeros(n+1)
      @simd for j in range(1,k)
          @simd for i in range(1,n+1)
              @inbounds temparray3[i,j] = AFArray(zeros(3))
          end
      end

      Bodycoordstemp = Bodycoords( AFArray(temparray2), AFArray(tempvector), temparray3 )
end
######################################
# MassData
function func_init_MassData(T::Type)
      tempfloat::T = 0.0     # Kappaleen tiheys

      MassDatatemp = MassData( tempfloat, tempfloat, tempfloat )
end
######################################
# Inertias
function func_init_Inertias(T::Type)
      temparray = zeros(T, 3, 3)

      Inertiastemp = Inertias( temparray, deepcopy(temparray), deepcopy(temparray), deepcopy(temparray) )
end
function func_init_Inertias(T::Type{StaticArray})
      temparray = @MMatrix(zeros(3, 3))

      Inertiastemp = Inertias( temparray, deepcopy(temparray), deepcopy(temparray), deepcopy(temparray) )
end
function func_init_Inertias(T::Type{AFArray})
      temparray = zeros(3, 3)

      Inertiastemp = Inertias( AFArray(temparray), AFArray(temparray), AFArray(temparray), AFArray(temparray) )
end
######################################
# StateVariables
function func_init_StateVariables(T::Type)
      tempvector3 = zeros(T, 3)
      tempvector4 = zeros(T, 4)

      StateVariablestemp = StateVariables( tempvector3, tempvector4, deepcopy(tempvector3), deepcopy(tempvector3) )
end
function func_init_StateVariables(T::Type{StaticArray})
      tempvector3 = @MVector(zeros(3))
      tempvector4 = @MVector(zeros(4))

      StateVariablestemp = StateVariables( tempvector3, tempvector4, deepcopy(tempvector3), deepcopy(tempvector3) )
end
function func_init_StateVariables(T::Type{AFArray})
      tempvector3 = zeros(3)
      tempvector4 = zeros(4)

      StateVariablestemp = StateVariables( AFArray(tempvector3), AFArray(tempvector4), AFArray(tempvector3), AFArray(tempvector3) )
end
######################################
# Derivaatat
function func_init_Derivaatat(T::Type)
      tempvector3 = zeros(T, 3)
      tempvector4 = zeros(T, 4)

      Derivaatattemp = Derivaatat( tempvector3, deepcopy(tempvector3), tempvector4, deepcopy(tempvector4) )
end
function func_init_Derivaatat(T::Type{StaticArray})
      tempvector3 = @MVector(zeros(3))
      tempvector4 = @MVector(zeros(4))

      Derivaatattemp = Derivaatat( tempvector3, deepcopy(tempvector3), tempvector4, deepcopy(tempvector4) )
end
function func_init_Derivaatat(T::Type{AFArray})
      tempvector3 = zeros(3)
      tempvector4 = zeros(4)

      Derivaatattemp = Derivaatat( AFArray(tempvector3), AFArray(tempvector3), AFArray(tempvector4), AFArray(tempvector4) )
end
######################################
# ApuArvot
function func_init_ApuArvot(T::Type)
      tempvector3 = zeros(T, 3)
      temparray = zeros(T, 3,3)
      temparray2 = zeros(T, 4,4)
      temparray3 = zeros(T, 3,4)

      ApuArvottemp = ApuArvot( temparray, deepcopy(temparray), tempvector3, deepcopy(tempvector3), deepcopy(temparray), temparray2, deepcopy(temparray2), temparray3 )
end
function func_init_ApuArvot(T::Type{StaticArray})
      tempvector3 = @MVector(zeros(3))
      temparray = @MMatrix(zeros(3,3))
      temparray2 = @MMatrix(zeros(4,4))
      temparray3 = @MMatrix(zeros(3,4))

      ApuArvottemp = ApuArvot( temparray, deepcopy(temparray), tempvector3, deepcopy(tempvector3), deepcopy(temparray), temparray2, deepcopy(temparray2), temparray3 )
end
function func_init_ApuArvot(T::Type{AFArray})
      tempvector3 = zeros(3)
      temparray = zeros(3,3)
      temparray2 = zeros(4,4)
      temparray3 = zeros(3,4)

      ApuArvottemp = ApuArvot( AFArray(temparray), AFArray(temparray), AFArray(tempvector3), AFArray(tempvector3), AFArray(temparray), AFArray(temparray2), AFArray(temparray2), AFArray(temparray3) )
end
######################################
# Voimat
function func_init_Voimat(T::Type)
      tempvector3 = zeros(T, 3)

      Voimattemp = Voimat( tempvector3, deepcopy(tempvector3), deepcopy(tempvector3), deepcopy(tempvector3) )
end
function func_init_Voimat(T::Type{StaticArray})
      tempvector3 = @MVector(zeros(3))

      Voimattemp = Voimat( tempvector3, deepcopy(tempvector3), deepcopy(tempvector3), deepcopy(tempvector3) )
end
function func_init_Voimat(T::Type{AFArray})
      tempvector3 = zeros(3)

      Voimattemp = Voimat( AFArray(tempvector3), AFArray(tempvector3), AFArray(tempvector3), AFArray(tempvector3) )
end
######################################
# KontaktiData
function func_init_KontaktiData(n::Integer, k::Integer, T::Type)
      tempfloat::T = 0.0
      tempvector3 = zeros(T, 3)
      temparray = zeros(T, n+1,k)
      temparray3 = Array(Array{T},n+1,k)
      @simd for j in range(1,k)
          @simd for i in range(1,n+1)
              @inbounds temparray3[i,j] = deepcopy(tempvector3)
          end
      end

      KontaktiDatatemp = KontaktiData( tempfloat, tempfloat, temparray3, deepcopy(temparray3), deepcopy(temparray3), temparray, deepcopy(temparray), deepcopy(temparray) )
end
function func_init_KontaktiData(n::Integer, k::Integer, T::Type{StaticArray})
      tempfloat::T = 0.0
      tempvector3 = @MVector(zeros(3))
      temparray = @MMatrix(zeros(n+1,k))
      temparray3 = Array(StaticArrays.MVector{3,Float64},n+1,k)
      @simd for j in range(1,k)
          @simd for i in range(1,n+1)
              @inbounds temparray3[i,j] = deepcopy(tempvector3)
          end
      end

      KontaktiDatatemp = KontaktiData( tempfloat, tempfloat, temparray3, deepcopy(temparray3), deepcopy(temparray3), temparray, deepcopy(temparray), deepcopy(temparray) )
end
function func_init_KontaktiData(n::Integer, k::Integer, T::Type{AFArray})
      tempfloat = 0.0
      temparray = zeros(n+1,k)
      temparray3 = Array(ArrayFire.AFArray{Float64},n+1,k)
      a = similar(temparray3)
      b = similar(temparray3)
      @simd for j in range(1,k)
          @simd for i in range(1,n+1)
              @inbounds temparray3[i,j] = AFArray(zeros(3))
              @inbounds a[i,j] = AFArray(zeros(3))
              @inbounds b[i,j] = AFArray(zeros(3))
          end
      end

      KontaktiDatatemp = KontaktiData( tempfloat, tempfloat, temparray3, a, b, AFArray(temparray), AFArray(temparray), AFArray(temparray) )
end
######################################
# BoundingBox
function func_init_BoundingBox(T::Type)
      tempvector2 = zeros(T, 2)
      temparray = zeros(T, 2, 3)

      BoundingBoxtemp = BoundingBox( tempvector2, deepcopy(tempvector2), deepcopy(tempvector2), temparray )
end
function func_init_BoundingBox(T::Type{StaticArray})
      tempvector2 = @MVector(zeros(2))
      temparray = @MMatrix(zeros(2, 3))

      BoundingBoxtemp = BoundingBox( tempvector2, deepcopy(tempvector2), deepcopy(tempvector2), temparray )
end
function func_init_BoundingBox(T::Type{AFArray})
      tempvector2 = zeros(2)
      temparray = zeros(2, 3)

      BoundingBoxtemp = BoundingBox( AFArray(tempvector2), AFArray(tempvector2), AFArray(tempvector2), AFArray(temparray) )
end
######################################
# KitkaData
function func_init_KitkaData(T::Type, T2::Type)
      tempfloat::T = 0.0
      tempint::T2 = 0

      KitkaDatatemp = KitkaData( tempfloat, tempfloat, tempfloat, tempfloat, tempfloat, tempfloat, tempint, tempfloat, tempfloat )
end
#################################################################################
# Init funktiot kappaleille
# Luodaan kappale, jolla on n,k pistettä ja pituus w
function func_init_body_empty(n::Integer, k::Integer, w::AbstractFloat, T::Type)
      pd = func_init_pd(n::Integer, k::Integer, w::AbstractFloat)
      body = func_init_Bodycoords(n::Integer, k::Integer, T::Type)
      world = func_init_Worldcoords(n::Integer, k::Integer, T::Type)
      md = func_init_MassData(typeof(w))
      J = func_init_Inertias(T::Type)
      sv = func_init_StateVariables(T::Type)
      dv = func_init_Derivaatat(T::Type)
      aux = func_init_ApuArvot(T::Type)
      f = func_init_Voimat(T::Type)
      knt = func_init_KontaktiData(n::Integer, k::Integer, T::Type)
      kd = func_init_KitkaData(typeof(w), typeof(n))
      bb = func_init_BoundingBox(T::Type)

      body_empty = RigidBody( pd, body, world, md, J, sv, dv, aux, f, knt, kd, bb )
      return body_empty
end
#################################################################################
# Kappaleiden alustusfunktiot
# Funktio, joka täyttää olemassa olevan kappaleen datasolut
function fill_rigidbody_pölli!(kpl::RigidBody)
      # Täyttää sylinterikoordinaatit pöllin muodolla R_alku halkaisijalla
      func_sektorit_serial!(kpl.body.R, kpl.body.Θ, kpl.pd.n, kpl.pd.k, R_alku)
      # Päivitetään lokaali karteesinen koordinaatisto vastaamaan sylinterikoordinaattien muotoa
      func_XYZ_pölli_rigidbody!(kpl.body.XYZ, kpl.body.R, kpl.body.Θ, kpl.pd.n, kpl.pd.k, kpl.pd.w)
      # päivitetään tiheys
      kpl.md.ρ = ρ
      # Lasketaan tilauus ja massa
      kpl.md.volume = pi*kpl.body.R[1,1]^2*kpl.pd.w
      kpl.md.massa = kpl.md.ρ*kpl.md.volume
      # Hitausmomentit
      func_MoI_cylinder!(kpl.J.body, kpl.md.massa, kpl.body.R[1,1], kpl.pd.w)
      kpl.J.body_inv = inv(kpl.J.body)
      # Eulerin parametrit ja rotaatiomatriisit. Oletuksena on että alussa lokaali koordinaatisto on samansuuntainen globaalin kanssa
      kpl.sv.q[1] = 1.0
      func_rotation_matrix!(kpl.aux.Rot_mat,kpl.sv.q)
      kpl.aux.Rot_mat_inv = inv(kpl.aux.Rot_mat)

      return nothing
end
