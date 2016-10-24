# Functions for modeling contact between bodies #
#################################################
type Contact
      A::Kappale # Kappale joka sisältää verteksin
      B::Kappale # Kappale joka sisältää pinnan

      p::AbstractArray # Kontaktipisteen koordinaatit globaalissa koordinaatistossa
      n::AbstractArray # Pinnan normaali, joka osoittaa ulospäin
      eA::AbstractArray # Reunan suuntavektori kappaleelle A
      eB::AbstractArray # Reunan suuntavektori kappaleelle B

      vf::Bool # true jos kontakti on verteksi/pinta
end

# funktio, joka laskee rangaistusvoiman
function func_F_penalty(kp, painuma, kv, v_relative, ki, Iterm, Nkontaktit)
      # kp = jousivakio
      # painuma = kontaktipisteiden välinen penetraatio
      # kv = Vaimennusvakio
      # v_relative = kontaktipisteiden välinen suhteellinen nopeus
      # ki = integraatiokerroinvakio
      # Iterm = Integraatiotermi, joka laskee maksimipainuman virheen summautumista simulaation aikana
      # Nkontaktit = kontaktipisteiden kokonaismäärä, joihin kohdistuu rangaistusvoima
      F_penalty = max( (kp*painuma - kv*v_relative + ki*Iterm) / Nkontaktit, 0.0)
end

# Funktio joka laskee kontaktipisteiden suhteellisen nopeuden
function func_v_relative(ṗa::AbstractArray, ṗb::AbstractArray, normal::AbstractArray)
      # kontaktinormaali osoittaa kappaleen a suuntaan
      v_relative = dot(normal, (ṗa - ṗb))
end

# Funktio joka laskee kontaktipisteen nopeuden
function func_p_dot(a::Kappale, point::AbstractArray)
      pdot = a.dv.ẋ + cross(a.aux.ω_world, (point - a.sv.x) )
end

# Funktio, joka laskee integraatiotermin
function func_Iterm(t_step::Integer, max_dyi::Array)
      λ::Float64 = 0.85
      Iterm::Float64 = 0.0

      @simd for i = collect(1:t_step-1)
            @inbounds Iterm += λ^(t_step-i-1) * max_dyi[i]
      end
      return Iterm
end

# Function that checks for contact between ground and bodies' vertices
function func_ground_contact!(Nbodies, contact_list)
      empty!(contact_list)

      for ii = collect(2:length(Nbodies))
            knt_pituus = length(contact_list)

            @simd for j in range(1,Nbodies[ii].pd.k)
                  @simd for i in range(1,Nbodies[ii].pd.n)
                        @inbounds begin
                              if (Nbodies[ii].world.X[i,j] < Nbodies[1].world.X[1] || Nbodies[ii].world.X[i,j] > Nbodies[1].world.X[3] || Nbodies[ii].world.Z[i,j]  < Nbodies[1].world.Z[1] || Nbodies[ii].world.Z[i,j]  > Nbodies[1].world.Z[2])
                                    #println("Ei oo sisällä!")
                                    #kontakti[1] = 0
                              else
                                    if Nbodies[ii].world.Y[i,j] < Nbodies[1].world.Y[1]
                                          #println("On sisällä!")
                                          #kontakti[1] = 1
                                          contact_point = Nbodies[ii].world.XYZ[i,j]
                                          contact_instance = Contact(Nbodies[ii], Nbodies[1], contact_point, [0.0; 1.0; 0.0], zeros(3), zeros(3), true)
                                          push!(contact_list, contact_instance)
                                    else
                                          #println("Ei oo sisällä!")
                                          #kontakti[1] = 0
                                    end
                              end
                        end
                  end
            end
            if length(contact_list) == knt_pituus
                  empty!(Nbodies[ii].knt.max_dyi)
            end
      end
      return nothing
end

# Function that handles ground contacts
function func_ground_collisions!(contact_list, Nbodies)
      ϵ = 0.65
      tolerance = 0.00001
      Nimpulse = 0
      Npenalty = 0
      impulse_contact = Array(Int,0)
      penalty_contact = Array(Int,0)
      F_penalty = zeros(3)

      # Go through contact list and divide it between collisions and resting contacts
      for i = collect(1:length(contact_list))
            ṗa = func_p_dot(contact_list[i].A, contact_list[i].p)
            ṗb = func_p_dot(contact_list[i].B, contact_list[i].p)
            vrel = func_v_relative(ṗa, ṗb, contact_list[i].n)
            if vrel < -tolerance #Impulse
                  push!(impulse_contact, i)
                  Nimpulse += 1
            elseif -tolerance <= vrel <= tolerance #Resting contact
                  push!(penalty_contact, i)
                  Npenalty += 1

                  push!(contact_list[i].A.knt.max_dyi, 0.0)
                  penetration = contact_list[i].B.world.Y[1,1] - contact_list[i].p[2]

                  if abs(penetration) > abs(contact_list[i].A.knt.max_dyi[end])
                        contact_list[i].A.knt.max_dyi[end] = penetration
                  end
            end
      end
      # Calculate the impulse contacts first
      if Nimpulse > 0
            for i = collect(1:Nimpulse)
                  ṗa = func_p_dot(contact_list[impulse_contact[i]].A, contact_list[impulse_contact[i]].p)
                  ṗb = func_p_dot(contact_list[impulse_contact[i]].B, contact_list[impulse_contact[i]].p)
                  vrel = func_v_relative(ṗa, ṗb, contact_list[impulse_contact[i]].n)

                  ra = contact_list[impulse_contact[i]].p - contact_list[impulse_contact[i]].A.sv.x
                  rb = contact_list[impulse_contact[i]].p - contact_list[impulse_contact[i]].B.sv.x

                  numerator = -(1.0 + ϵ)*vrel
                  # Denominator
                  # body A
                  term1 = cross(ra, contact_list[impulse_contact[i]].n)
                  term2 = contact_list[impulse_contact[i]].A.J.world_inv*term1
                  term3 = cross(term2,ra)
                  term4 = dot(contact_list[impulse_contact[i]].n,term3)
                  # body B
                  term1B = cross(rb, contact_list[impulse_contact[i]].n)
                  term2B = contact_list[impulse_contact[i]].B.J.world_inv*term1B
                  term3B = cross(term2B,rb)
                  term4B = dot(contact_list[impulse_contact[i]].n,term3B)

                  j_impulse = numerator / ( contact_list[impulse_contact[i]].A.md.massa_inv + contact_list[impulse_contact[i]].B.md.massa_inv + term4 + term4B)

                  # Impulse force
                  F_impulse = contact_list[impulse_contact[i]].n * j_impulse / Nimpulse
                  # Apply impulse to bodies
                  contact_list[impulse_contact[i]].A.sv.P += F_impulse
                  contact_list[impulse_contact[i]].B.sv.P -= F_impulse

                  contact_list[impulse_contact[i]].A.sv.L += cross(ra, F_impulse)
                  contact_list[impulse_contact[i]].B.sv.L -= cross(rb, F_impulse)
            end
            # Update velocities
            for i  in length(Nbodies)
                  Nbodies[i].dv.ẋ = Nbodies[i].sv.P*Nbodies[i].md.massa_inv
                  Nbodies[i].aux.ω_world = Nbodies[i].J.world_inv*Nbodies[i].sv.L
            end
      end
      # Calculate penaltymethods
      for i = collect(1:Npenalty)
            ṗa = func_p_dot(contact_list[penalty_contact[i]].A, contact_list[penalty_contact[i]].p)
            ṗb = func_p_dot(contact_list[penalty_contact[i]].B, contact_list[penalty_contact[i]].p)
            vrel = func_v_relative(ṗa, ṗb, contact_list[penalty_contact[i]].n)

            ra = contact_list[penalty_contact[i]].p - contact_list[penalty_contact[i]].A.sv.x
            rb = contact_list[penalty_contact[i]].p - contact_list[penalty_contact[i]].B.sv.x
            penetration = contact_list[penalty_contact[i]].B.world.Y[1,1] - contact_list[penalty_contact[i]].p[2]

            t_step = length(contact_list[penalty_contact[i]].A.knt.max_dyi)
            Iterm = func_Iterm(t_step, contact_list[penalty_contact[i]].A.knt.max_dyi)
            # Calculate the penalty force coefficients
            kp = 180.0*( contact_list[penalty_contact[i]].A.md.massa + contact_list[penalty_contact[i]].B.md.massa )
            kv = 90.0*( contact_list[penalty_contact[i]].A.md.massa + contact_list[penalty_contact[i]].B.md.massa )
            ki = 4.0*( contact_list[penalty_contact[i]].A.md.massa + contact_list[penalty_contact[i]].B.md.massa )
            # Calculate the penalty force
            F_penalty[2] = func_F_penalty(kp, penetration, kv, vrel, ki, Iterm, Npenalty)
            # Apply penalty force and torque to bodies
            contact_list[penalty_contact[i]].A.f.force_world += F_penalty
            contact_list[penalty_contact[i]].B.f.force_world -= F_penalty

            contact_list[penalty_contact[i]].A.f.torque_world += cross(ra, F_penalty)
            contact_list[penalty_contact[i]].B.f.torque_world -= cross(rb, F_penalty)
      end
      return nothing
end
