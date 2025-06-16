for datuminfo in DATUMINFO
  wgsdatum = Symbol("WGS84", datuminfo)
  refdatum = Symbol("REF", datuminfo)
  funcdatum = datuminfo

 @eval $funcdatum(::Type{<:T}=Float64) where T<:AbstractFloat = convert(T, $refdatum[])
 @eval $funcdatum(::Type{<:Real})  = convert(Float64, $refdatum[])
end


# setdatum! function to set the datum parameters
@eval function setdatum!(datum::Datum=WGS84Latest)::Nothing where Datum

     _EARTH = ellipsoid(datum)

    _MAJORAXIS = majoraxis(_EARTH) |> er-> uconvert(km,er) |> ustrip
    _MINORAXIS = minoraxis(_EARTH) |> er-> uconvert(km,er) |> ustrip
    _ECCENTRICITY² = eccentricity²(_EARTH) |> ustrip
    _ECCENTRICITY = sqrt(_ECCENTRICITY²)   |> ustrip
    _NORMALIZEMINORAXIS = _MINORAXIS / _MAJORAXIS
    _COMPLECCENTRICITY² = 1 - _ECCENTRICITY²
    $(SETDATUM...)
    return nothing
end

# setmajoraxis! function to set the major axis
@eval function setmajoraxis!(majoraxis::T=REFMAJORAXIS[])::Nothing where T<:AbstractFloat
      _MAJORAXIS = majoraxis
      _MINORAXIS = REFMINORAXIS[] # use the current minor axis
      if _MAJORAXIS < _MINORAXIS
        @warn("Major axis must be greater than or equal to minor axis. But got: $majoraxis < $(REFMINORAXIS[]). Keeping current values ($REFMAJORAXIS, $REFMINORAXIS).")
        return nothing
      end
      _ECCENTRICITY² = 1 - (_MINORAXIS / _MAJORAXIS)^2
      _ECCENTRICITY = sqrt(_ECCENTRICITY²)
      _NORMALIZEMINORAXIS = _MINORAXIS / _MAJORAXIS
      _COMPLECCENTRICITY² = 1 - _ECCENTRICITY²
      $(SETDATUM...)
    return nothing
end

# setminoraxis! function to set the minor axis
@eval function setminoraxis!(minoraxis::T=REFMINORAXIS[])::Nothing where T<:AbstractFloat
      _MINORAXIS = minoraxis
      _MAJORAXIS = REFMAJORAXIS[] # use the current major axis
      if _MINORAXIS > _MAJORAXIS
        @warn("Minor axis must be less than or equal to major axis. But got: $minoraxis > $(REFMAJORAXIS[]). Keeping current values ($REFMAJORAXIS, $REFMINORAXIS).")
        return nothing
      end
      _ECCENTRICITY² = 1 - (_MINORAXIS / _MAJORAXIS)^2
      _ECCENTRICITY = sqrt(_ECCENTRICITY²)
      _NORMALIZEMINORAXIS = _MINORAXIS / _MAJORAXIS
      _COMPLECCENTRICITY² = 1 - _ECCENTRICITY²
      $(SETDATUM...)
    return nothing
end

# seteccentricity²! function to set the eccentricity squared
@eval function seteccentricity²!(eccentricity²::T=REFECCENTRICITY²[])::Nothing where T<:AbstractFloat
      _ECCENTRICITY² = eccentricity²
      if _ECCENTRICITY² < 0 || _ECCENTRICITY² > 1
        @warn("Eccentricity squared must be in the range [0, 1]. But got: $eccentricity². Keeping current values ($(REFECCENTRICITY²[])).")
        return nothing
      end
      _ECCENTRICITY = sqrt(_ECCENTRICITY²)
      _MAJORAXIS = REFMAJORAXIS[] # use the current major axis
      _MINORAXIS = _MAJORAXIS * sqrt(1 - _ECCENTRICITY²)
      _NORMALIZEMINORAXIS = _MINORAXIS / _MAJORAXIS
      _COMPLECCENTRICITY² = 1 - _ECCENTRICITY²
      $(SETDATUM...)
    return nothing
end


@inline function clampangle(θ::T, θmin::T, θmax::T)::T where T<:AbstractFloat
    if θmin < θmax
      (θmin<θ<θmax) && return θ
    else
      (θmin<θ || θ<θmax) && return θ
    end
    distθmin = mod(θmin-θ,360)
    distθmax = mod(θ-θmax,360)
    return distθmin < distθmax ? θmin : θmax
end


@inline function intersectionrayray(pointx1::T,pointy1::T, directionx1::T, directiony1::T,
                                    pointx2::T,pointy2::T, directionx2::T, directiony2::T) where T<:AbstractFloat

    Δx = pointx2 - pointx1
    Δy = pointy2 - pointy1

    det = -directionx1*directiony2 + directiony1*directionx2

    det= abs(det) < eps(T) ? eps(T)*sign(det) : det

    t = (Δy*directionx2 - Δx*directiony2)/det
    h = (Δy*directionx1 - Δx*directiony1)/det

    return t,h
end


@inline function rotation_matrix(angle::Real,nx::Real,ny::Real)
  cosθ = cosd(angle)
  sinθ = sind(angle)
  return cosθ*nx+sinθ*ny, -sinθ*nx+cosθ*ny
end

@inline function nadir_angle_normal(nx::Real,ny::Real,angle::Real)
  norm_ = hypot(nx, ny)
  nx,ny = nx/norm_, ny/norm_

  rotation_matrix(angle, nx, ny)
end

@inline function limb_angle_normal(nx::Real,ny::Real,angle::Real)
  nadir_angle_normal(ny, -nx, angle)
end

"""
    infologger(filename::String)
Create a logger that writes log messages to a file.
This logger uses the `FormatLogger` to format messages and write them to the specified file.
"""
function infologger(filename::String)
  FormatLogger(filename;) do io, args
    println(io, args.message)
  end
end
