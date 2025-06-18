using Makie
using Makie.MakieCore
using Yagami
# bring in missing Makie methods required for block definition
using Makie: make_block_docstring
const Rect2d = Rect2{Float64}
Makie.@Block RayAxis <: Makie.AbstractAxis begin
  scene::Scene
  targetlimits::Observable{Rect2d}
  finallimits::Observable{Rect2d}
  mouseeventhandle::Makie.MouseEventHandle
  scrollevents::Observable{Makie.ScrollEvent}
  keysevents::Observable{Makie.KeysEvent}
  interactions::Dict{Symbol, Tuple{Bool, Any}}
  elements::Dict{Symbol, Any}
  transform_func::Observable{Any}
  inv_transform_func::Observable{Any}

  @attributes begin
    xscale =identity
    yscale =identity
    # layout observables for Block
    "The horizontal alignment of the block in its suggested bounding box."
    halign = :center
    "The vertical alignment of the block in its suggested bounding box."
    valign = :center
    "The width setting of the block."
    width = Makie.Auto()
    "The height setting of the block."
    height = Makie.Auto()
    "Control if the parent layout can adjust to this block's width."
    tellwidth::Bool = true
    "The align mode of the block in its parent GridLayout."
    alignmode = Makie.Inside()

    # datum
    "Major axis of the erarh, default is $(Yagami.YagamiCore.WGS84MAJORAXIS) km."
    majoraxis = Yagami.YagamiCore.WGS84MAJORAXIS
    "Minor axis of the earth, default is $(Yagami.YagamiCore.WGS84MAJORAXIS) km."
    minoraxis = Yagami.YagamiCore.WGS84MAJORAXIS
    # levels
    "Level set for the ray tracing. It is a vector of heights in km."
    levels = Yagami.RayTracing.create_hlevelset(Float32;hmin= 4,hmax= 120,levels=10)
    # radii
    "Radii for the ray tracing. It is a vector of angles in degrees."
    radii  = Yagami.RayTracing.create_radii(Float32;θmin=0,θmax=350,radii=36)
    "Visualization of the levels on a logarithmic scale. Default is `false`."
    simpleviz::Bool =


    "Controls if the y axis goes upwards (false) or downwards (true)"
    yreversed::Bool = false
    "Controls if the x axis goes rightwards (false) or leftwards (true)"
    xreversed::Bool = false
    "The relative margins added to the autolimits in x direction."
    xautolimitmargin::Tuple{Float64,Float64} = (0.05f0, 0.05f0)
    "The relative margins added to the autolimits in y direction."
    yautolimitmargin::Tuple{Float64,Float64} = (0.05f0, 0.05f0)
    "The limits that the user has manually set. They are reinstated when calling `reset_limits!` and are set to nothing by `autolimits!`. Can be either a tuple (xlow, xhigh, ylow, high) or a tuple (nothing_or_xlims, nothing_or_ylims). Are set by `xlims!`, `ylims!` and `limits!`."
    limits = (nothing, nothing)
    "The forced aspect ratio of the axis. `nothing` leaves the axis unconstrained, `DataAspect()` forces the same ratio as the ratio in data limits between x and y axis, `AxisAspect(ratio)` sets a manual ratio."
    aspect = Makie.DataAspect()
    autolimitaspect = nothing

    # appearance controls
    "The set of fonts which text in the axis should use.s"
    fonts = (; regular = "TeX Gyre Heros Makie")
    "The axis title string."
    title = ""
    "The font family of the title."
    titlefont = :bold
    "The title's font size."
    titlesize::Float64 = @inherit(:fontsize, 16f0)
    "The gap between axis and title."
    titlegap::Float64 = 4f0
    "Controls if the title is visible."
    titlevisible::Bool = true
    "The horizontal alignment of the title."
    titlealign::Symbol = :center
    "The color of the title"
    titlecolor::RGBAf = @inherit(:textcolor, :black)
    "The axis title line height multiplier."
    titlelineheight::Float64 = 1
    "The axis subtitle string."
    subtitle = ""
    "The font family of the subtitle."
    subtitlefont = :regular
    "The subtitle's font size."
    subtitlesize::Float64 = @inherit(:fontsize, 16f0)
    "The gap between subtitle and title."
    subtitlegap::Float64 = 0
    "Controls if the subtitle is visible."
    subtitlevisible::Bool = true
    "The color of the subtitle"
    subtitlecolor::RGBAf = @inherit(:textcolor, :black)
    "The axis subtitle line height multiplier."
    subtitlelineheight::Float64 = 1


    "The xlabel string."
    xlabel = ""
    "The ylabel string."
    ylabel = ""
    "The font family of the xlabel."
    xlabelfont = :regular
    "The font family of the ylabel."
    ylabelfont = :regular
    "The color of the xlabel."
    xlabelcolor::RGBAf = @inherit(:textcolor, :black)
    "The color of the ylabel."
    ylabelcolor::RGBAf = @inherit(:textcolor, :black)
    "The font size of the xlabel."
    xlabelsize::Float64 = @inherit(:fontsize, 16f0)
    "The font size of the ylabel."
    ylabelsize::Float64 = @inherit(:fontsize, 16f0)
    "Controls if the xlabel is visible."
    xlabelvisible::Bool = true
    "Controls if the ylabel is visible."
    ylabelvisible::Bool = true
    "The padding between the xlabel and the ticks or axis."
    xlabelpadding::Float64 = 3f0
    "The padding between the ylabel and the ticks or axis."
    ylabelpadding::Float64 = 5f0 # xlabels usually have some more visual padding because of ascenders, which are larger than the hadvance gaps of ylabels
    "The xlabel rotation in radians."
    xlabelrotation = Makie.automatic
    "The ylabel rotation in radians."
    ylabelrotation = Makie.automatic

    "The x (longitude) ticks - can be a vector or a Makie tick finding algorithm."
    xticks = Makie.automatic
    "The y (latitude) ticks - can be a vector or a Makie tick finding algorithm."
    yticks = Makie.automatic

    "Format for x (longitude) ticks."
    xtickformat = Makie.automatic
    "Format for y (latitude) ticks."
    ytickformat = Makie.automatic
    "The font family of the xticklabels."
    xticklabelfont = :regular
    "The font family of the yticklabels."
    yticklabelfont = :regular
    "The color of xticklabels."
    xticklabelcolor::RGBAf = @inherit(:textcolor, :black)
    "The color of yticklabels."
    yticklabelcolor::RGBAf = @inherit(:textcolor, :black)
    "The font size of the xticklabels."
    xticklabelsize::Float64 = @inherit(:fontsize, 16f0)
    "The font size of the yticklabels."
    yticklabelsize::Float64 = @inherit(:fontsize, 16f0)
    "Controls if the xticklabels are visible."
    xticklabelsvisible::Bool = true
    "Controls if the yticklabels are visible."
    yticklabelsvisible::Bool = true
    "The space reserved for the xticklabels."
    xticklabelspace::Union{Makie.Automatic, Float64} = Makie.automatic
    "The space reserved for the yticklabels."
    yticklabelspace::Union{Makie.Automatic, Float64} = Makie.automatic
    "The space between xticks and xticklabels."
    xticklabelpad::Float64 = 5f0
    "The space between yticks and yticklabels."
    yticklabelpad::Float64 = 5f0
    "The counterclockwise rotation of the xticklabels in radians."
    xticklabelrotation::Float64 = 0f0
    "The counterclockwise rotation of the yticklabels in radians."
    yticklabelrotation::Float64 = 0f0
    "The horizontal and vertical alignment of the xticklabels."
    xticklabelalign::Union{Makie.Automatic, Tuple{Symbol, Symbol}} = Makie.automatic
    "The horizontal and vertical alignment of the yticklabels."
    yticklabelalign::Union{Makie.Automatic, Tuple{Symbol, Symbol}} = Makie.automatic
    "The size of the xtick marks."
    xticksize::Float64 = 6f0
    "The size of the ytick marks."
    yticksize::Float64 = 6f0
    "Controls if the xtick marks are visible."
    xticksvisible::Bool = true
    "Controls if the ytick marks are visible."
    yticksvisible::Bool = true
    "The alignment of the xtick marks relative to the axis spine (0 = out, 1 = in)."
    xtickalign::Float64 = 0f0
    "The alignment of the ytick marks relative to the axis spine (0 = out, 1 = in)."
    ytickalign::Float64 = 0f0
    "The width of the xtick marks."
    xtickwidth::Float64 = 1f0
    "The width of the ytick marks."
    ytickwidth::Float64 = 1f0
    "The color of the xtick marks."
    xtickcolor::RGBAf = RGBf(0, 0, 0)
    "The color of the ytick marks."
    ytickcolor::RGBAf = RGBf(0, 0, 0)
    # "The width of the axis spines."
    # spinewidth::Float64 = 1f0
    "Controls if the x grid lines are visible."
    xgridvisible::Bool = true
    "Controls if the y grid lines are visible."
    ygridvisible::Bool = true
    "The width of the x grid lines."
    xgridwidth::Float64 = 1f0
    "The width of the y grid lines."
    ygridwidth::Float64 = 1f0
    "The color of the x grid lines."
    xgridcolor::RGBAf = RGBAf(0, 0, 0, 0.5)
    "The color of the y grid lines."
    ygridcolor::RGBAf = RGBAf(0.0, 0, 0, 0.5)
    "The linestyle of the x grid lines."
    xgridstyle = nothing
    "The linestyle of the y grid lines."
    ygridstyle = nothing
    "Controls if minor ticks on the x axis are visible"
    xminorticksvisible::Bool = false
    "The alignment of x minor ticks on the axis spine"
    xminortickalign::Float64 = 0f0
    "The tick size of x minor ticks"
    xminorticksize::Float64 = 4f0
    "The tick width of x minor ticks"
    xminortickwidth::Float64 = 1f0
    "The tick color of x minor ticks"
    xminortickcolor::RGBAf = :black
    "The tick locator for the x minor ticks"
    xminorticks = IntervalsBetween(2)
    "Controls if minor ticks on the y axis are visible"
    yminorticksvisible::Bool = false
    "The alignment of y minor ticks on the axis spine"
    yminortickalign::Float64 = 0f0
    "The tick size of y minor ticks"
    yminorticksize::Float64 = 4f0
    "The tick width of y minor ticks"
    yminortickwidth::Float64 = 1f0
    "The tick color of y minor ticks"
    yminortickcolor::RGBAf = :black
    "The tick locator for the y minor ticks"
    yminorticks = IntervalsBetween(2)
    "Controls if the x minor grid lines are visible."
    xminorgridvisible::Bool = false
    "Controls if the y minor grid lines are visible."
    yminorgridvisible::Bool = false
    "The width of the x minor grid lines."
    xminorgridwidth::Float64 = 1f0
    "The width of the y minor grid lines."
    yminorgridwidth::Float64 = 1f0
    "The color of the x minor grid lines."
    xminorgridcolor::RGBAf = RGBAf(0, 0, 0, 0.05)
    "The color of the y minor grid lines."
    yminorgridcolor::RGBAf = RGBAf(0, 0, 0, 0.05)
    "The linestyle of the x minor grid lines."
    xminorgridstyle = nothing
    "The linestyle of the y minor grid lines."
    yminorgridstyle = nothing
    # "Controls if the axis spine is visible."
    # spinevisible::Bool = true
    # "The color of the axis spine."
    # spinecolor::RGBAf = :black
    # spinetype::Symbol = :geospine
    "The button for panning."
    panbutton::Makie.Mouse.Button = Makie.Mouse.right
    "The key for limiting panning to the x direction."
    xpankey::Makie.Keyboard.Button = Makie.Keyboard.x
    "The key for limiting panning to the y direction."
    ypankey::Makie.Keyboard.Button = Makie.Keyboard.y
    "The key for limiting zooming to the x direction."
    xzoomkey::Makie.Keyboard.Button = Makie.Keyboard.x
    "The key for limiting zooming to the y direction."
    yzoomkey::Makie.Keyboard.Button = Makie.Keyboard.y

    "Locks interactive panning in the x direction."
    xpanlock::Bool = false
    "Locks interactive panning in the y direction."
    ypanlock::Bool = false
    "Locks interactive zooming in the x direction."
    xzoomlock::Bool = false
    "Locks interactive zooming in the y direction."
    yzoomlock::Bool = false
    "Controls if rectangle zooming affects the x dimension."
    xrectzoom::Bool = true
    "Controls if rectangle zooming affects the y dimension."
    yrectzoom::Bool = true

    xaxisposition::Symbol = :bottom
    yaxisposition::Symbol = :left
  end
end



@recipe(AltitudeTrace) do scene
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

@recipe(RayTrace) do scene
  Attributes(
    :logscale => false
    :majoraxis => MAJORAXIS(),
    :minoraxis => MINORAXIS(),
    :levels => Base.LogRange(4,120,10),
    :radii  => LinRange(0,350,36),
  )
end

@inline function __ellipse(h,θ,majoraxis,minoraxis)
  @inline N₀(θ) = majoraxis/sqrt(1-minoraxis^2*sind(θ)^2)
  b=minoraxis/majoraxis
  b² = b^2
  Point2f((N₀(θ)+h)*cosd(θ), (N₀(θ)*b²+h)*sind(θ))
end


function Makie.plot!(
  sc::RayTrace{<:Tuple{<:AbstractVector{<:Real}, <:AbstractVector{<:Real}}})

  X= sc[1]
  Y= sc[2]

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



"""
    raytracing_grid!(sc::Axis, angles::AbstractVector{<:Real}, altitudes::AbstractVector{<:Real},
                     majoraxis::Real = WGS84MAJORAXIS,
                     minoraxis::Real = WGS84MINORAXIS,color=:black,linewidth=3,kwargs...)
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
