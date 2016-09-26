# Datatypes used for the simulation
# Rigid body datatype
type RigidBody
     massa::Float64
     # Hitausmomentit
     J_body::Array{Float64}
     J_body_inv::Array{Float64}
     J_world::Array{Float64}
     J_world_inv::Array{Float64}

     # State variables
     x::Vector{Float64}       # World coordinates of center of mass
     q::Vector{Float64}       # Quaternions
     P::Vector{Float64}       # Linear momentum
     L::Vector{Float64}       # Angular momentum

     # Auxiliary variables
     Rot_mat::Array{Float64}  # Rotation matrix
     v::Vector{Float64}       # Linear velocity of center of mass in world coordinates

     # Computed quantities
     force::Vector{Float64}   # Forces affecting the body in world coordinates
     torque::Vector{Float64}  # Torques acting on the body in world coordinates
end
