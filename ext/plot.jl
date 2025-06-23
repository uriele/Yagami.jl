@recipe(AltitudeTrace) do scene,args...
  Attributes(
    :tqmarker => :star6,
    :tqcolor  => :green,
    :tqsize  => 10,
  )
end

function Makie.plot!(
  sc::AltitudeTrace{<:Tuple{<:AbstractVector{<:Real}, <:AbstractVector{<:Real}}})

  angles= sc[1]
  altitudes = sc[2]

  tangentquote = Observable(Point2f)
  linesegs     =Observable(Point2f[])

  function update_plot!(angles, altitudes)

    empty!(linesegs[])

    _,idx=findmin(altitudes)

    tangentquote[] = Point2f(angles[idx], altitudes[idx])
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
