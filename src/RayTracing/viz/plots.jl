using Makie
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

"""
    `raytracing_grid!(ax, angles, altitudes;majoraxis = WGS84MAJORAXIS, minoraxis = WGS84MINORAXIS,color=:black,linewidth=3,kwargs...)`
Draws a grid of rays on an axis, where each ray is defined by a constant angle and altitude.
The rays are drawn as vertical lines at constant angles and horizontal lines at constant altitudes.
The grid is based on an ellipsoid defined by the major and minor axes.
"""
function raytracing_grid!(sc::Axis, angles::AbstractVector{<:Real}, altitudes::AbstractVector{<:Real},
                         majoraxis::Real = WGS84MAJORAXIS,
                         minoraxis::Real = WGS84MINORAXIS,color=:black,linewidth=3,kwargs...)

    a = majoraxis
    b = minoraxis
    e² = 1 - (b / a)^2  # no `[]` needed unless `a`, `b` are Observables

    θmin, θmax = extrema(angles)
    hmin, hmax = extrema(altitudes)

    @inline N₀(θ) = a / sqrt(1 - e² * sind(θ)^2)

    @inline ellipsepoints(θ, h) = Point2f(
        (N₀(θ) + h) * cosd(θ),
        (N₀(θ) * (1 - e²) + h) * sind(θ)
    )

    # Draw vertical lines at constant angle between min/max altitude
    @inbounds for θ in angles
        lines!(sc, [ellipsepoints(θ, hmin), ellipsepoints(θ, hmax)];
               color=color, linewidth=linewidth, kwargs...)
    end

    # Draw horizontal lines at constant altitude between θmin and θmax
    @inbounds for h in altitudes
        thetas = range(θmin, θmax, length=200)
        points = [ellipsepoints(θ, h) for θ in thetas]
        lines!(sc, points;
               color=color, linewidth=linewidth, kwargs...)
    end

    return sc
end
