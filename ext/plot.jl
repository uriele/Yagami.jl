using Makie

@recipe(RayTrace) do scene
  Attributes(
    :tqmarker => :star6,
    :tqcolor  => :green,
    :tqsize  => 10,
  )
end

function Makie.plot!(
  sc::RayTrace{<:Tuple{<:AbstractVector{<:Real}, <:AbstractVector{<:Real}}})

  angles= sc[1]
  altitudes = sc[2]

  tangentquote = Observable(Point2f)
  linesegs     =Observable(Point2f[])

  function update_plot!(angles, altitudes)

    empty!(linesegs[])

    _,idx=findmin(altitudes)

    tangentquote = Point2f(angles[idx], altitudes[idx])
    for (angle,altitude) in zip(angles, altitudes)
      push!(linesegs[], Point2f(angle, altitude))
    end
  end
    Makie.Observables.onany(update_plot!, angles, altitudes)

    update_plot!(angles[], altitudes[])

    lines!(sc,linesegs)
    scatter!(sc, tangentquote; marker=sc.tqmarker,color=sc.tqcolor, markersize=sc.tqsize)

    sc
end


figure = Figure()


angles = LinRange(0,180,100)
altitudes = @. 100 - 50*sind(angles)
altitudes2 = @. 80 - 50*sind(angles)
angles
raytrace(collect(angles), altitudes;tqcolor=:red,tqsize=20)
raytrace!(collect(angles), altitudes2)
