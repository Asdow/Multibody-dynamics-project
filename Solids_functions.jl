# Functions to build solids of type RigidBody
#######################################
# Creates a "cylinder" that is divided into n,k points, w is length of cylinder
# Function that populates a body's datastructs for a cylinder
function fill_rigidbody_cylinder!(kpl::RigidBody, Radius::Float64, ρ::Float64)
      # Täyttää sylinterikoordinaatit pöllin muodolla R_alku halkaisijalla
      func_sektorit_cylinder!(kpl.body.R, kpl.body.Θ, kpl.pd.n, kpl.pd.k, R_alku)
      # Päivitetään lokaali karteesinen koordinaatisto vastaamaan sylinterikoordinaattien muotoa
      func_XYZ_cylinder_rigidbody!(kpl.body.XYZ, kpl.body.R, kpl.body.Θ, kpl.pd.n, kpl.pd.k, kpl.pd.w)
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

# Creates a cube that is h high, w wide and l deep
function fill_rigidbody_cube!(body::RigidBody, h::Float64, w::Float64, l::Float64, ρ::Float64)
      # Päivitetään lokaali karteesinen koordinaatisto vastaamaan kuution muotoa
      func_XYZ_cube_rigidbody!(body.body.XYZ, h, w, l, body.pd.n, body.pd.k)
      #func_z_midpoint_rigidbody!(body.XYZ_body, n, k, w)
      # päivitetään tiheys
      body.md.ρ = ρ
      # Lasketaan tilauus ja massa
      body.md.volume = func_vol_cube(h, w, l)
      body.md.massa = body.md.ρ*body.md.volume
      body.md.massa_inv = inv(body.md.massa)
      # Hitausmomentit
      func_MoI_cube!(body.J.body, body.md.massa, h, w, l)
      body.J.body_inv = inv(body.J.body)
      # Eulerin parametrit ja rotaatiomatriisit. Oletuksena on että alussa lokaali koordinaatisto on samansuuntainen globaalin kanssa
      body.sv.q[1] = 1.0
      func_rotation_matrix!(body.aux.Rot_mat,body.sv.q)
      body.aux.Rot_mat_inv = inv(body.aux.Rot_mat)

      return nothing
end

# Creates a Ground plane that is w wide and l deep and at height h from global origin
function fill_rigidbody_Ground!(kpl::RigidBody, w, l, h)
      kpl.pd.w = w
#=      kpl.world.X = [-w/2.0; -w/2.0; w/2.0; w/2.0; -w/2.0]
      kpl.world.Y = [h; h; h; h; h]
      kpl.world.Z = [-l/2.0; l/2.0; l/2.0; -l/2.0; -l/2.0]
      @simd for i in range(1,5)
            kpl.world.XYZ[i] = [kpl.world.X[i]; kpl.world.Y[i]; kpl.world.Z[i]]
      end
=#
      kpl.sv.x[2] = h
      body_X = [-w/2.0; -w/2.0; w/2.0; w/2.0; -w/2.0]
      body_Y = [0.0; 0.0; 0.0; 0.0; 0.0]
      body_Z = [-l/2.0; l/2.0; l/2.0; -l/2.0; -l/2.0]
      @simd for i in range(1,5)
            kpl.body.XYZ[i] = [body_X[i]; body_Y[i]; body_Z[i]]
      end
      # Eulerin parametrit ja rotaatiomatriisit. Oletuksena on että alussa lokaali koordinaatisto on samansuuntainen globaalin kanssa
      kpl.sv.q[1] = 1.0
      func_rotation_matrix!(kpl.aux.Rot_mat,kpl.sv.q)
      kpl.aux.Rot_mat_inv = inv(kpl.aux.Rot_mat)
      func_body_to_world(kpl)

      # Set ground's coefficients for the penalty method
      kpl.knt.k = 5.0e3
      kpl.knt.c = 1.5e3
      kpl.knt.ki = 1.0e2
      
      return nothing
end
