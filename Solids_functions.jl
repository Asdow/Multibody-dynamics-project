# Functions to build solids of type RigidBody 
#######################################
# Creates a "cylinder" that is divided into n,k points, w is length of cylinder
function fill_rigidbody_cylinder!(body::RigidBody, n::Int64, k::Int64, w::Float64, Radius::Float64, ρ::Float64)
      # Täyttää sylinterikoordinaatit pöllin muodolla Radius halkaisijalla
      func_sektorit_cylinder!(body.R, body.theta, n, k, Radius)
      # Päivitetään lokaali karteesinen koordinaatisto vastaamaan sylinterikoordinaattien muotoa
      func_XYZ_cylinder_rigidbody!(body.XYZ_body, body.R, body.theta, n, k, w)
      # päivitetään tiheys
      body.rho = ρ
      # Lasketaan tilauus ja massa
      body.volume = func_vol_cylinder(Radius, w)
      body.massa = body.rho*body.volume
      # Hitausmomentit
      func_MoI_cylinder!(body.J_body, body.massa, Radius, w)
      body.J_body_inv = inv(body.J_body)
      # Eulerin parametrit ja rotaatiomatriisit. Oletuksena on että alussa lokaali koordinaatisto on samansuuntainen globaalin kanssa
      body.q[1] = 1.0
      func_rotation_matrix!(body.Rot_mat,body.q)
      body.Rot_mat_inv = inv(body.Rot_mat)

      return nothing
end

# Creates a cube that is h high, w wide and l deep
function fill_rigidbody_cube!(body::RigidBody, n::Int64, k::Int64, h::Float64, w::Float64, l::Float64, ρ::Float64)
      # Päivitetään lokaali karteesinen koordinaatisto vastaamaan kuution muotoa
      func_XYZ_cube_rigidbody!(body.XYZ_body, h, w, l, n, k)
      #func_z_midpoint_rigidbody!(body.XYZ_body, n, k, w)
      # päivitetään tiheys
      body.rho = ρ
      # Lasketaan tilauus ja massa
      body.volume = func_vol_cube(h, w, l)
      body.massa = body.rho*body.volume
      # Hitausmomentit
      func_MoI_cube!(body.J_body, body.massa, h, w, l)
      body.J_body_inv = inv(body.J_body)
      # Eulerin parametrit ja rotaatiomatriisit. Oletuksena on että alussa lokaali koordinaatisto on samansuuntainen globaalin kanssa
      body.q[1] = 1.0
      func_rotation_matrix!(body.Rot_mat,body.q)
      body.Rot_mat_inv = inv(body.Rot_mat)

      return nothing
end
