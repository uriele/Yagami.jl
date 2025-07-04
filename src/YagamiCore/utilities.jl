for datuminfo in DATUMINFO
  wgsdatum = Symbol("WGS84", datuminfo)
  refdatum = Symbol("REF", datuminfo)
  funcdatum = datuminfo

 @eval $funcdatum(::Type{<:T}=Float64) where T<:AbstractFloat = convert(T, $refdatum[])
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

@inline difference_angles(θ1, θ2) = mod(θ1 - θ2 + 180, 360) - 180
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


@inline function rotation_matrix(angle::T,nx::T,ny::T) where T
  cosθ = cosd(angle)
  sinθ = sind(angle)
  return cosθ*nx+sinθ*ny, -sinθ*nx+cosθ*ny
end

@inline function nadir_angle_normal(nx::T,ny::T,angle::T) where T
  norm_ = hypot(nx, ny)
  nx,ny = nx/norm_, ny/norm_

  rotation_matrix(angle, nx, ny)
end

@inline function limb_angle_normal(nx::T,ny::T,angle::T) where T
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

@inline numshort(n::T;digits=3,padding=10) where T<:Real = lpad(@sprintf("%.*e",digits,n),padding+digits)
@inline numshortf(n::T;digits=3,padding=10) where T<:Real = lpad(@sprintf("%.*f",digits,n),padding+digits)
@inline numshort(n::Int;digits=3,padding=10) = lpad(@sprintf("%*i",digits,n),padding+digits)
@inline textshort(n::String;digits=3,padding=10) = string(lpad(n,padding+digits))

@doc """
    `file_to_array(filename,header,type=Float64,spacer=" ")`

Read a file and convert its contents into a vector of type `type`. The file is expected to have a header of `header` lines, which will be skipped.
The contents of the file are split by the specified `spacer`, and each element is parsed to the specified `type`.

# Arguments:
- `filename::String`: The path to the file to be read.
- `header::Int=0`: The number of header lines to skip in the file (default is 0).
- `type::Type{T}=Float64`: The type to which the elements of the file will be parsed (default is `Float64`).
- `spacer::S=" "`: The string used to split the contents of the file (default is a single space).

# Returns:
- `Vector{T}`: An array of type `T` containing the parsed elements from the file.

# Example:
```julia
julia> file_to_array("data.txt", header=1, type=Float64, spacer=",")
[1.0, 2.0, 3.0, 4.0]
```
"""
@inline function file_to_array(filename::String,header::Int=0,type::Type{T}=Float64,spacer::S=" ") where {T,S<:AbstractString}
  open(filename,"r") do io
    lines=readlines(io)
    lines=lines[header+1:end] # skip header lines
    data = join(lines, spacer) |> d-> split(d, spacer) |>
    d-> filter(x -> !isempty(x),d) |>
    d-> map(x -> parse(T, x),d)
  end
end


@inline function vrm_to_namedtuple(file,header,type::Type{T}=Float64,spacer=" ") where T<:AbstractFloat  #find the number of pollutants using regex
  open(file,"r") do io
    lines=readlines(io)[header+1:end] # skip header lines
    regex=Regex("^\\s{3,}[1-9][0-9]?\\s+([a-zA-Z0-9]+)")
    nlines=length(lines)

    idx= Int[]
    pollutant=String[]

    for (i,line) in enumerate(lines)
      if occursin(regex,line)
        m=match(regex,line)
        push!(idx,i)
        push!(pollutant,m.captures[1])
      end
    end

    @inline _convert_to_array(lines,spacer) = join(lines, spacer) |>
    d-> split(d, spacer) |>
    d-> filter(x -> !isempty(x),d) |>
    d-> map(x -> parse(T, x),d)

    vmrs=Dict{String,Vector{Float64}}(pollutant[end] =>
    _convert_to_array(lines[idx[end]+1:end],spacer))
    for i in eachindex(idx[1:end-1])
      # get the lines between two indices and
      vmrs[pollutant[i]] = _convert_to_array(lines[idx[i]+1:idx[i+1]-1],spacer)
    end
    (; (Symbol(k) => v for (k, v) in vrms)...)
  end
end
