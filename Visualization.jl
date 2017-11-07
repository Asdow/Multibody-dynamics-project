function init_screen()
    window = gl.glscreen(); # Avataan ikkuna johon plotataan.
    @async gl.renderloop(window); # Renderöidään kuvaaja. async, koska renderloop blokkaa
    window
end

# funktio jolla testataan visualisointi yms.
function gltest(pölli)
    window = init_screen();
    poss = fposs(pölli);
    # Luodaan Array{GLAbstraction.Context{GLAbstraction.DeviceUnit}}, jotta voidaan tallentaa pöllin visualisointi yhteen muuttujaan
    pöllivis = GLAbstraction.Context{GLAbstraction.DeviceUnit}[];
    # Piirretään pöllin pisteet pieninä palloina
    spheres = gl.visualize((gt.Sphere{Float32}(gt.Point3f0(0.0), 0.002f0), poss), boundingbox=nothing);
    # gl._view(spheres, window, camera=:perspective)
    # Luodaan yhteydet pisteiden välille viivojen piirtoa varten
    # Indeksit kertoo minkä pisteiden välille viivat piirretään.
    indeksit = indexointi(pölli);
    # Luodaan visualisointi
    lines = gl.visualize(poss, :linesegment, thickness=0.5f0, indices=indeksit, boundingbox=nothing);
    # gl._view(lines, window, camera=:perspective);
    # pusketaan visualisoinnit vektoriin
    push!(pöllivis, spheres)
    push!(pöllivis, lines)
    # Ladataan visualisointi, jotta se näkyy
    gl._view(pöllivis, window, camera=:perspective)
    return window, pöllivis
end
# w, vis = gltest(pölli);
"""
    fposs(kpl::Kappale)
Luodaan kappaleen pisteiden asemista GeometryTypes:n Points tyypit
"""
function fposs(kpl::Kappale)
    poss = gt.GeometryTypes.Point{3,Float32}[];
    for j in 1:size(kpl.coords.XYZb,2)
        for i in 1:size(kpl.coords.XYZb,1)
            push!(poss, gt.Point3f0(kpl.coords.XYZb[i,j]))
        end
    end
    return poss
end
"""
    indexointi_sylinteri(pölli::Kappale)
Luo indeksit kappaleen world pisteiden välille viivojen piirtoa varten
"""
function indexointi_sylinteri(kpl::Kappale)
    indeksit = Array{UInt32}(0);
    prof_pituus = size(kpl.world.XYZ,1)
    # Viivat profiilien pisteiden välillä
    for i in range(0,prof_pituus)
        for j in range(0,size(kpl.world.XYZ,2))
            push!(indeksit, Cuint(i+j*prof_pituus))
        end
    end
    # Viivat profiilin vastakkaiden pisteiden välillä (Vain pöllin päädyille)
    puolikas = prof_pituus÷2
    for n in (0,size(kpl.world.XYZ,2)-1)
        for i in range(0,puolikas)
            for j in range(0,2)
                push!(indeksit, Cuint(i+j*puolikas+n*prof_pituus))
            end
        end
    end

    return indeksit
end
"""
    indexointi_cube(pölli::Kappale)
Luo indeksit kappaleen world pisteiden välille viivojen piirtoa varten
"""
function indexointi_cube(kpl::Kappale)
      indeksit = Array{UInt32}(0);
      prof_pituus = size(kpl.coords.XYZ,1)
      # Viivat profiilien pisteiden välillä
      for i in range(0,prof_pituus)
            for j in range(0,size(kpl.coords.XYZ,2))
            push!(indeksit, Cuint(i+j*prof_pituus))
            end
      end
      # päädyt
      for j in range(0,size(kpl.coords.XYZ,2))
            for i in range(0,prof_pituus)
                  push!(indeksit, Cuint(i+j*prof_pituus))
            end
      end
      push!(indeksit, Cuint(1))
      push!(indeksit, Cuint(2))
      push!(indeksit, Cuint(3))
      push!(indeksit, Cuint(0))
      push!(indeksit, Cuint(5))
      push!(indeksit, Cuint(6))
      push!(indeksit, Cuint(7))
      push!(indeksit, Cuint(4))
      return indeksit
end
"""
    rotationmatrix(rotmax33::T) where {T<:StaticArrays.MArray{Tuple{3,3},Float64,2,9}}
Muuttaa kappaleen orientaation GLVisualizen vaatimaan muotoon.
"""
function rotationmatrix(rotmax33::T) where {T<:StaticArrays.MArray{Tuple{3,3},Float64,2,9}}
    T0 = zero(Float32)
    T1 = one(Float32)
    gla.Mat4f0(
        rotmax33[1,1],  rotmax33[2,1],  rotmax33[3,1],  T0,
        rotmax33[1,2],  rotmax33[2,2],  rotmax33[3,2],  T0,
        rotmax33[1,3],  rotmax33[2,3],  rotmax33[3,3],  T0,
        T0, T0, T0, T1,
    )
end
"""
    translationmatrix(t::StaticArrays.MVector{3, T}) where {T<:AbstractFloat}
Muuttaa kappaleen aseman GLVisualizen vaatimaan muotoon.
"""
function translationmatrix(t::sa.MVector{3, T}) where {T<:AbstractFloat}
    T0, T1 = zero(Float32), one(Float32)
    gla.Mat4f0(
        T1,  T0,  T0,  T0,
        T0,  T1,  T0,  T0,
        T0,  T0,  T1,  T0,
        t[1],t[2],t[3],T1,
    )
end
"""
    transformation(kpl::Kappale)
Muuttaa kappaleen aseman ja orientaation GLVisualizen vaatimaan muotoon.
"""
function transformation(kpl::Kappale)
    T0, T1 = zero(Float32), one(Float32)
    rotmax33 = kpl.aux.R;
    t = kpl.sv.x;
    gla.Mat4f0(
    rotmax33[1,1],  rotmax33[2,1],  rotmax33[3,1],  T0,
    rotmax33[1,2],  rotmax33[2,2],  rotmax33[3,2],  T0,
    rotmax33[1,3],  rotmax33[2,3],  rotmax33[3,3],  T0,
    t[1], t[2], t[3], T1,
    )
end
"""
    body_vis(body::Kappale)
Creates a body's visualization. Returns a GLAbstraction.Context
"""
function body_vis(body::Kappale)
      # Array that holds the body's visualizations
      bodyvis = gla.Context{gla.DeviceUnit}[];
      push!(bodyvis, gl.visualize(body.sh.mesh, :lines, thickness = 1f0, color = col.RGBA(1f0, 0f0, 0f0, 0.8f0), boundingbox=nothing))
      # push!(bodyvis, gl.visualize(body.sh.mesh, boundingbox=nothing))
      # body reference coordinate axes
      push!(bodyvis, axes(0.8f0, 1.0f0))
      # transform all visualizations to body's global location and orientation
      gl.set_arg!(bodyvis, :model, transformation(body))
      return bodyvis
end
"""
    body_vis!(bodyvis::Array{GLAbstraction.Context{GLAbstraction.DeviceUnit},1}, body::T) where {T<:Kappale}
Updates a body's visualization inplace.
"""
function body_vis!(bodyvis::Array{GLAbstraction.Context{GLAbstraction.DeviceUnit},1}, body::T) where {T<:Kappale}
    gl.set_arg!(bodyvis, :model, transformation(body))
    return nothing
end
"""
    body_vis!(bodyvis::Array{Array{GLAbstraction.Context{GLAbstraction.DeviceUnit},1},1}, bodies::Array{T,1}) where {T<:Kappale}
Updates bodies' visualizations inplace.
"""
function body_vis!(bodyvis::Array{Array{GLAbstraction.Context{GLAbstraction.DeviceUnit},1},1}, bodies::Array{T,1}) where {T<:Kappale}
    for i in 1:length(bodies)
        body_vis!(bodyvis[i], bodies[i])
    end
    return nothing
end
"""
    origin()
Draws the global coordinate axes.
"""
function origin()
      o = axes(1.0f0, 2.0f0);
      return o::GLAbstraction.Context{GLAbstraction.DeviceUnit}
end
"""
    axes(len::Float32, t::Float32)
Creates differently colored linesegments in the direction of 3 axes. Len adjusts the length and t the thickness.
"""
function axes(len::Float32, t::Float32)
      # Points describing the line segment positions.
      ax = [gt.Point3f0(0), gt.Point3f0(len,0,0), gt.Point3f0(0), gt.Point3f0(0,len,0), gt.Point3f0(0), gt.Point3f0(0,0,len)];
      # Colors for the line segments. length must be equal to positions.
      cols = [col.RGBA(1f0,0f0,0f0,1f0), col.RGBA(1f0,0f0,0f0,1f0), col.RGBA(0f0,1f0,0f0,1f0), col.RGBA(0f0,1f0,0f0,1f0), col.RGBA(0f0,0f0,1f0,1f0), col.RGBA(0f0,0f0,1f0,1f0)];
      axesvis = gl.visualize(ax, :linesegment, thickness=t, color=cols, boundingbox=nothing);
      return axesvis::GLAbstraction.Context{GLAbstraction.DeviceUnit}
end


# Function which sets an argument of a Context/RenderObject.
# If multiple RenderObjects are supplied, it'll try to set the same argument in all
# of them.
# gl.set_arg!(robj::Context, sym::Symbol, value)
# Esim
# gl.set_arg!(lines, :model, transformation(pölli))

# Way to get the attributedata/keywords from a visualization.
# eg. bodyvis.children[].uniforms
