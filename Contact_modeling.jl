# Functions for modeling contact between bodies #
#################################################
abstract type ForceElement{T<:Real} end
struct ForceE2{T<:Real} <: ForceElement{T}
    A::RigidBody{T}
    B::RigidBody{T}
    F::sa.SVector{3,T}
    T::sa.SVector{3,T}
    Tb::sa.SVector{3,T}
end
struct ForceE1{T<:Real} <: ForceElement{T}
    A::RigidBody{T}
    F::sa.SVector{3,T}
    T::sa.SVector{3,T}
end
# Contacts are modeled through penalty method based on:
# Evan Drumwright.
# "A Fast and Stable Penalty Method for Rigid Body Simulation."
# IEEE Trans. on Visualization and Computer Graphics, vol. 14, no. 1, pp. 231-240, January/February, 2008.
# funktio, joka laskee rangaistusvoiman
function F_penalty(kp, painuma, kv, v_relative, ki, Iterm, Nkontaktit)
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
function vreln(ṗa, ṗb, normal)
      # Contact normal points to body A.
      v_relative = dot(normal, (ṗa - ṗb))
end
function vreln(ṗa, normal)
      # Contact normal points to body A.
      v_relative = dot(normal, ṗa)
end

# Funktio joka laskee kontaktipisteen nopeuden
function pdot(a::Kappale, point)
      pdot = a.dv.ẋ + cross(a.dv.ω, (point - a.sv.x) )
end
function pdot!(pdot, a::Kappale, point)
      pdot[:] = a.dv.ẋ + cross(a.dv.ω, (point - a.sv.x) )
end

# Funktio, joka laskee integraatiotermin
function Iterm(max_dyi)
      λ::Float64 = 0.85; # Forgetting factor
      Iterm::Float64 = 0.0;
      t_step = countnz(max_dyi)

      for i = 1:(t_step-1)
            @inbounds Iterm += λ^(t_step-i-1) * max_dyi[i]
      end
      return Iterm
end

function collision(a::Kappale, b::Kappale)
    # Get point of collision and penetration depth
    point, depth = getcollposdepth(b,p)
    # Assume plane is always stationary, so vrel = ṗa
    ṡ = pdot(b, point)
    ṡn = vreln(ṡ, p.normal)
    It = Iterm(b.knt.max_dyi)

    Fp = F_penalty(b.knt.coeff.k, depth, b.knt.coeff.c, ṡn, b.knt.coeff.ki, It, 1)
    Fpvec = Fp.*p.normal # Penaltyforce points towards sphere
    # Calculate moment affecting the body, as Collshape isn't necessarily on ref origin.
    ra = point - b.sv.x;
    T = cross(ra, Fpvec)
    # Create a Force struct that only points to body b.
    ForceE1(b, Fpvec, T)
end

function collision(b::Kappale, p::Plane)
    # Get point of collision and penetration depth
    point, depth = getcollposdepth(b,p)
    # Assume plane is always stationary, so vrel = ṗa
    ṡ = pdot(b, point)
    ṡn = vreln(ṡ, p.normal)
    It = Iterm(b.knt.max_dyi)

    Fp = F_penalty(b.knt.coeff.k, depth, b.knt.coeff.c, ṡn, b.knt.coeff.ki, It, 1)
    Fpvec = Fp.*p.normal # Penaltyforce points towards sphere
    # Calculate moment affecting the body, as Collshape isn't necessarily on ref origin.
    ra = point - b.sv.x;
    T = cross(ra, Fpvec)
    # Create a Force struct that only points to body b.
    ForceE1(b, Fpvec, T)
end


function resolvecoll(Rsys, list, p::Plane)
    Flist = ForceElement{eltype(Rsys)}[];
    for i in list
        push!(Flist, collision(Rsys.bodies[i], p))
    end
    return Flist
end
function resolvecoll(Rsys, Flist)
    Flist = ForceElement{eltype(Rsys)}[];
    for i in list
        push!(Flist, collision(Rsys.bodies[i[1]], Rsys.bodies[i[2]]))
    end
    return Flist
end
