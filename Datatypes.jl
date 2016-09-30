# Datatypes used for the simulation
###################################### #
# Rigid body datatype
type RigidBody
      # Kappaleen koordinaattidata
      n::Int64              # Kappaleen pisteiden lukumäärä
      k::Int64              # Kappaleen pisteiden lukumäärä
      R::Array{Float64}       # Kappaleen muoto lokaalissa sylinterikoordinaatistossa, säde
      Θ::Vector{Float64}   # Kappaleen muoto lokaalissa sylinterikoordinaatistossa, kulma
      XYZ_body::Array{Array}  # Kappaleen pisteiden asema lokaalissa kooridnaatistossa
      XYZ_world::Array{Array}  # Kappaleen pisteiden asema globaalissa kooridnaatistossa
      # Kappaleen massa ja tilavuusdata
      rho::Float64      # Kappaleen tiheys
      volume::Float64   # Kappaleen tilavuus
      massa::Float64    # Kappaleen massa
      # Hitausmomentit
      J_body::Array{Float64}
      J_body_inv::Array{Float64}
      J_world::Array{Float64}
      J_world_inv::Array{Float64}
      # State variables
      x::Vector{Float64}        # Massakeskipisteen asema globaalissa koordinaatistossa
      q::Vector{Float64}        # Eulerin parametrit
      P::Vector{Float64}        # Lineaarinen liikemäärä globaalissa koordinaatistossa
      L::Vector{Float64}        # Angulaarinen liikemäärä
      # Apuarvot
      Rot_mat::Array{Float64}         # Rotaatiomatriisi
      Rot_mat_inv::Array{Float64}     # Rotaatiomatriisin käänteismatriisi
      a::Vector{Float64}              # Massakeskipisteen lineaarinen kiihtyvyys globaalissa koordinaatistossa
      v::Vector{Float64}              # Massakeskipisteen lineaarinen nopeus globaalissa koordinaatistossa
      omega_body::Vector{Float64}     # Kappaleen kulmanopeus lokaalissa koordinaatistossa
      omega_world::Vector{Float64}    # Kappaleen kulmanopeus globaalissa koordinaatistossa
      qa::Vector{Float64}             # Eulerin parametrien kiihtyvyys
      qdot::Vector{Float64}           # Eulerin parametrien muutosnopeus
      omega_tilde::Array{Float64}     # 3x3 skew symmetric omega_tilde matriisi
      omega_tilde4::Array{Float64}    # 4x4 skew symmetric omega_tilde matriisi
      E::Array{Float64}               # 4x4 E matriisi
      E_1::Array{Float64}             # 3x4 E matriisi
            # Lasketut arvot
      force::Vector{Float64}         # Kappaleeseen kohdistuvat voimat globaalissa koordinaatistossa
      torque_body::Vector{Float64}   # Kappaleeseen kohdistuvat väännöt lokaalissa koordinaatistossa
      torque_world::Vector{Float64}  # Kappaleeseen kohdistuvat väännöt globaalissa koordinaatistossa
end

# Init funktiot kappaleille
# Luodaan tyhjä kappale
function func_init_body_empty(n, k)
      tempfloat = 0.0
      tempvector3 = zeros(3)
      tempvector4 = zeros(4)
      tempvectorn = zeros(n+1)
      temparray = zeros(3,3)
      temparray2 = zeros(n+1,k)
      temparray3 = Array(Array{Float64},n+1,k)
      temparray4 = zeros(4,4)
      temparray5 = zeros(3,4)
      @simd for j in range(1,k)
          @simd for i in range(1,n+1)
              @inbounds temparray3[i,j] = deepcopy(tempvector3)
          end
      end

      body = RigidBody(n, k, deepcopy(temparray2), deepcopy(tempvectorn), deepcopy(temparray3), deepcopy(temparray3), tempfloat, tempfloat, tempfloat, deepcopy(temparray), deepcopy(temparray), deepcopy(temparray), deepcopy(temparray), deepcopy(tempvector3), deepcopy(tempvector4), deepcopy(tempvector3), deepcopy(tempvector3), deepcopy(temparray), deepcopy(temparray), deepcopy(tempvector3), deepcopy(tempvector3), deepcopy(tempvector3), deepcopy(tempvector3), deepcopy(tempvector4), deepcopy(tempvector4), deepcopy(temparray), deepcopy(temparray4), deepcopy(temparray4), deepcopy(temparray5), deepcopy(tempvector3), deepcopy(tempvector3), deepcopy(tempvector3) )
      return body
end
