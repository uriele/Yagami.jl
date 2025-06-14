
"""
    `Ray2D{T}`

  Data structure for ray tracing in 2D space.
# Arguments
  - `pointx::T`: x-coordinate of the ray's starting point in km. Default is 0.0.
  - `pointy::T`: y-coordinate of the ray's starting point in km. Default is 0.0.
  - `directionx::T`: x-component of the ray's normalized direction vector. Default is 0.0.
  - `directiony::T`: y-component of the ray's normalized direction vector. Default is 1.0.
  - `i::Int`: index of the wedge in which the ray is currently located. Default is 0.
  - `j::Int`: index of the sub-wedge in which the ray is currently located. Default is 0.

"""
struct Ray2D{T<:AbstractFloat}
    pointx::T
    pointy::T
    directionx::T
    directiony::T
    i:: Int
    j:: Int

    function Ray2D{T}(pointx::T=0.0, pointy::T=0.0, directionx::T=0.0, directiony::T=1.0, i::Int=0, j::Int=0) where T<:AbstractFloat
      dirnorm = hypot(directionx, directiony)
      directionx /= dirnorm
      directiony /= dirnorm
      new(pointx, pointy, directionx, directiony, i, j)
    end
end
"""
  $SIGNATURES

Data structures for ray tracing in 2D. `Ray2D` represents a ray in 2D space with a point and direction and the current wedge index

# Arguments
  - `pointx::T`: x-coordinate of the ray's starting point in km. Default is 0.0.
  - `pointy::T`: y-coordinate of the ray's starting point in km. Default is 0.0.
  - `directionx::T`: x-component of the ray's normalized direction vector. Default is 0.0.
  - `directiony::T`: y-component of the ray's normalized direction vector. Default is 1.0.
  - `i::Int`: index of the wedge in which the ray is currently located. Default is 0.
  - `j::Int`: index of the sub-wedge in which the ray is currently located. Default is 0.

"""
Ray2D(pointx::T=0.0, pointy::T=0.0, directionx::T=0.0, directiony::T=1.0, i::Int=0, j::Int=0) where T<:AbstractFloat =
 Ray2D{T}(pointx, pointy, directionx, directiony, i, j)
function Ray2D(pointx::Real,pointy::Real, directionx::Real, directiony::Real, i::Int=0, j::Int=0)
 args=promote(pointx, pointy, directionx, directiony)
 Ray2D{eltype(args)}(args..., i, j)
end

Base.show(io::IO, ::MIME"text/plain", r::Ray2D{T}) where {T} = begin
  println(io, "Rays2D{",T,"}($(__getpointx(r)),$(__getpointy(r)))km")
end
Base.show(io::IO,  r::Ray2D{T}) where {T} = begin
  println(io, "($(round(__getpointx(r),digits=2)),$(round(__getpointy(r),digits=2)))")
end

"""
 $SIGNATURES
Get the point and direction of a `Ray2D` as a 2D vector.
Returns a `SVector{2,T}` where `T` is the type of the ray's coordinates.
"""
getpoint(r::Ray2D{T}) where {T} = SVector{2,T}(r.pointx, r.pointy)
"""
  $SIGNATURES
Get the direction of a `Ray2D` as a 2D vector.
Returns a `SVector{2,T}` where `T` is the type of the ray's direction.
"""
getdirection(r::Ray2D{T}) where {T} = SVector{2,T}(r.directionx, r.directiony)

"""
  $SIGNATURES
Get the x or y coordinate of the point of a `Ray2D`.
- `sym::Symbol`: Either `:x` or `:y` to specify which coordinate to return.
"""
getpoint(r::Ray2D{T}, sym::Symbol) where {T} = if sym == :x
    r.pointx
elseif sym == :y
    r.pointy
end
"""
  $SIGNATURES
Get the x or y coordinate of the point of a `Ray2D`.
- `ind::Int`: Either `1` for x-coordinate or `2` for y-coordinate.
"""
getpoint(r::Ray2D{T}, ind::Int) where {T} = if ind == 1
    r.pointx
elseif ind == 2
    r.pointy
end
"""
  $SIGNATURES
Get the x or y component of the direction of a `Ray2D`.
- `sym::Symbol`: Either `:x` or `:y` to specify which component to return.
"""
getdirection(r::Ray2D{T}, sym::Symbol) where {T} = if sym == :x
    r.directionx
elseif sym == :y
    r.directiony
end
"""
  $SIGNATURES
Get the x or y component of the direction of a `Ray2D`.
- `ind::Int`: Either `1` for x-component or `2` for y-component.
"""
getdirection(r::Ray2D{T}, ind::Int) where {T} = if ind == 1
    r.directionx
elseif ind == 2
    r.directiony
end

"""
  $SIGNATURES

Get the i or j index of the wedge of a `Ray2D`.
- `sym::Symbol`: Either `:i` or `:j` to specify which index to return.
"""
getwedgeindex(r::Ray2D{T}, sym::Symbol) where {T} = if sym == :i
    r.i
elseif sym == :j
    r.j
end
"""
  $SIGNATURES
Get the i or j index of the wedge of a `Ray2D`.
- `ind::Int`: Either `1` for i-index or `2` for j-index.
"""
getwedgeindex(r::Ray2D{T}, ind::Int) where {T} = if ind == 1
    r.i
elseif ind == 2
    r.j
end


"""
`ResultRay{T}`

Data structure to hold the result of a ray intersection in 2D space.
This structure contains the intersection point, direction, index of the wedge,
length of the intersection, altitude, and azimuth of the ray.
# Arguments
  - `intersection_pointx::T`: x-coordinate of the intersection point in km.
  - `intersection_pointy::T`: y-coordinate of the intersection point in km.
  - `intersection_directionx::T`: x-component of the ray's normalized direction vector at the intersection.
  - `intersection_directiony::T`: y-component of the ray's normalized direction vector at the intersection.
  - `intersection_index_i::Int`: index of the wedge in which the ray intersects.
  - `intersection_index_j::Int`: index of the sub-wedge in which the ray intersects.
  - `intersection_length::T`: length of the intersection in km.
  - `intersection_altitude::T`: altitude at the intersection point in km.
  - `intersection_azimuth::T`: azimuth angle at the intersection point in radians.

"""
struct ResultRay{T}
    intersection_pointx::T
    intersection_pointy::T
    intersection_directionx::T
    intersection_directiony::T
    intersection_index_i::Int
    intersection_index_j::Int
    intersection_length::   T
    intersection_altitude:: T
    intersection_azimuth::  T

    ResultRay{T}(directionx::T=0.0, directiony::T=0.0, pointx::T=0.0, pointy::T=0.0,
                  index_i::Int=0, index_j::Int=0, length::T=0.0,
                  altitude::T=0.0, azimuth::T=0.0) where T =  new(pointx, pointy, directionx, directiony, index_i, index_j, length, altitude, azimuth)
end


"""
  $SIGNATURES
Get the intersection point and direction of a `ResultRay` as a 2D vector.
Returns a `SVector{2,T}` where `T` is the type of the ray's coordinates.
"""
getpoint(r::ResultRay{T}) where {T} = SVector{2,T}(r.intersection_pointx, r.intersection_pointy)
"""
  $SIGNATURES
Get the intersection direction of a `ResultRay` as a 2D vector.
Returns a `SVector{2,T}` where `T` is the type of the ray's direction.
"""
getdirection(r::ResultRay{T}) where {T} = SVector{2,T}(r.intersection_directionx, r.intersection_directiony)


"""
  $SIGNATURES
Get the x or y coordinate of the intersection point of a `ResultRay`.
- `sym::Symbol`: Either `:x` or `:y` to specify which coordinate to return.
"""
getpoint(r::ResultRay{T}, sym::Symbol) where {T} = if sym == :x
    r.intersection_pointx
elseif sym == :y
    r.intersection_pointy
end
"""
  $SIGNATURES
Get the x or y coordinate of the intersection point of a `ResultRay`.
- `ind::Int`: Either `1` for x-coordinate or `2` for y-coordinate.
"""
getpoint(r::ResultRay{T}, ind::Int) where {T} = if ind == 1
    r.intersection_pointx
elseif ind == 2
    r.intersection_pointy
end


"""
  $SIGNATURES
Get the x or y component of the intersection direction of a `ResultRay`.
- `sym::Symbol`: Either `:x` or `:y` to specify which component to return.
"""
getdirection(r::ResultRay{T}, sym::Symbol) where {T} = if sym == :x
    r.intersection_directionx
elseif sym == :y
    r.intersection_directiony
end
"""
  $SIGNATURES
Get the x or y component of the intersection direction of a `ResultRay`.
- `ind::Int`: Either `1` for x-component or `2` for y-component.
"""
getdirection(r::ResultRay{T}, ind::Int) where {T} = if ind == 1
    r.intersection_directionx
elseif ind == 2
    r.intersection_directiony
end
Base.show(io::IO, r::ResultRay{T}) where {T} = begin
  println(io, "ResultRay{",T,"}($(r.intersection_pointx),$(r.intersection_pointy))")
end

"""
  $SIGNATURES
Get the i or j index of the wedge of a `ResultRay`.
- `sym::Symbol`: Either `:i` or `:j` to specify which index to return.
"""
getwedgeindex(r::ResultRay{T}, sym::Symbol) where {T} = if sym == :i
    r.intersection_index_i
elseif sym == :j
    r.intersection_index_j
end
"""
  $SIGNATURES
Get the i or j index of the wedge of a `ResultRay`.
- `ind::Int`: Either `1` for i-index or `2` for j-index.
"""
getwedgeindex(r::ResultRay{T}, ind::Int) where {T} = if ind == 1
    r.intersection_index_i
elseif ind == 2
    r.intersection_index_j
end


"""
  $SIGNATURES
Get the length, altitude, and azimuth of a `ResultRay`. Returns the length of the intersection in km.
"""
getlength(r::ResultRay{T}) where {T} = r.intersection_length
""" $SIGNATURES
Get the altitude and azimuth of a `ResultRay`. Returns the altitude in km.
"""
getaltitude(r::ResultRay{T}) where {T} = r.intersection_altitude
"""
  $SIGNATURES
Get the azimuth of a `ResultRay`. Returns the azimuth angle in degrees.
"""
getazimuth(r::ResultRay{T}) where {T} = r.intersection_azimuth


############################################################################################
# Inline getters and setters for Ray2D and ResultRay (internal use)
#############################################################################################

@inline __getpointx(r::Ray2D{T}) where {T} = r.pointx
@inline __getpointy(r::Ray2D{T}) where {T} = r.pointy
@inline __getdirectionx(r::Ray2D{T}) where {T} = r.directionx
@inline __getdirectiony(r::Ray2D{T}) where {T} = r.directiony
@inline __geti(r::Ray2D{T}) where {T} = r.i
@inline __getj(r::Ray2D{T}) where {T} = r.j



const Ray2Ds{T} = StructArray{Ray2D{T}, 1}
@inline __getpointx(r::Ray2Ds{T}, idx::Int) where {T} = r.pointx[idx]
@inline __getpointy(r::Ray2Ds{T}, idx::Int) where {T} = r.pointy[idx]
@inline __getdirectionx(r::Ray2Ds{T}, idx::Int) where {T} = r.directionx[idx]
@inline __getdirectiony(r::Ray2Ds{T}, idx::Int) where {T} = r.directiony[idx]
@inline __geti(r::Ray2Ds{T}, idx::Int) where {T} = r.i[idx]
@inline __getj(r::Ray2Ds{T}, idx::Int) where {T} = r.j[idx]

@inline __setpointx!(r::Ray2Ds{T}, idx::Int, x::T) where {T} = r.pointx[idx] = x
@inline __setpointy!(r::Ray2Ds{T}, idx::Int, y::T) where {T} = r.pointy[idx] = y
@inline __setdirectionx!(r::Ray2Ds{T}, idx::Int, x::T) where {T} = r.directionx[idx] = x
@inline __setdirectiony!(r::Ray2Ds{T}, idx::Int, y::T) where {T} = r.directiony[idx] = y
@inline __seti!(r::Ray2Ds{T}, idx::Int, i::Int) where {T} = r.i[idx] = i
@inline __setj!(r::Ray2Ds{T}, idx::Int, j::Int) where {T} = r.j[idx] = j


############################################################################################
# note: RayResult are immutable, so setters are not needed
############################################################################################
@inline __getpointx(r::ResultRay{T}) where {T} = r.intersection_pointx
@inline __getpointy(r::ResultRay{T}) where {T} = r.intersection_pointy
@inline __getdirectionx(r::ResultRay{T}) where {T} = r.intersection_directionx
@inline __getdirectiony(r::ResultRay{T}) where {T} = r.intersection_directiony
@inline __getlength(r::ResultRay{T}) where {T} = r.intersection_length
@inline __getaltitude(r::ResultRay{T}) where {T} = r.intersection_altitude
@inline __getazimuth(r::ResultRay{T}) where {T} = r.intersection_azimuth
@inline __geti(r::ResultRay{T}) where {T} = r.intersection_index_i
@inline __getj(r::ResultRay{T}) where {T} = r.intersection_index_j
##############

const ResultRays{T}= StructArray{ResultRay{T}, 1}

@inline __getpointx(r::ResultRays{T}     ,idx::Int) where {T} = r.intersection_pointx[idx]
@inline __getpointy(r::ResultRays{T}     ,idx::Int) where {T} = r.intersection_pointy[idx]
@inline __getdirectionx(r::ResultRays{T} ,idx::Int) where {T} = r.intersection_directionx[idx]
@inline __getdirectiony(r::ResultRays{T} ,idx::Int) where {T} = r.intersection_directiony[idx]
@inline __getlength(r::ResultRays{T}     ,idx::Int) where {T} = r.intersection_length[idx]
@inline __getaltitude(r::ResultRays{T}   ,idx::Int) where {T} = r.intersection_altitude[idx]
@inline __getazimuth(r::ResultRays{T}    ,idx::Int) where {T} = r.intersection_azimuth[idx]
@inline __geti(r::ResultRays{T}          ,idx::Int) where {T} = r.intersection_index_i[idx]
@inline __getj(r::ResultRays{T}          ,idx::Int) where {T} = r.intersection_index_j[idx]
##############

@inline __setpointx!(r::ResultRays{T}    ,idx::Int, x::T) where {T} = r.intersection_pointx[idx] = x
@inline __setpointy!(r::ResultRays{T}    ,idx::Int, y::T) where {T} = r.intersection_pointy[idx] = y
@inline __setdirectionx!(r::ResultRays{T},idx::Int, x::T) where {T} = r.intersection_directionx[idx] = x
@inline __setdirectiony!(r::ResultRays{T},idx::Int, y::T) where {T} = r.intersection_directiony[idx] = y
@inline __setlength!(r::ResultRays{T}    ,idx::Int, length::T) where {T} = r.intersection_length[idx] = length
@inline __setaltitude!(r::ResultRays{T}  ,idx::Int, altitude::T) where {T} = r.intersection_altitude[idx] = altitude
@inline __setazimuth!(r::ResultRays{T}   ,idx::Int, azimuth::T) where {T} = r.intersection_azimuth[idx] = azimuth
@inline __seti!(r::ResultRays{T}         ,idx::Int,i::Int) where {T} = r.intersection_index_i[idx] = i
@inline __setj!(r::ResultRays{T}         ,idx::Int,j::Int) where {T} = r.intersection_index_j[idx] = j
