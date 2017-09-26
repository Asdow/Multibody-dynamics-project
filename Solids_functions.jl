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


# Creates a Ground plane that is w wide and l deep and at height h from global origin
function create_plane(w, l, h)
      heightfield = zeros(Float32, w, l)
      fill!(heightfield, h)
      return heightfield
end
