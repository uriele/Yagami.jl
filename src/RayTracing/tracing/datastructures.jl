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
mutable struct Ray2D{T<:AbstractFloat}
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

@inline __getpointx(r::Ray2D{T}) where {T} = getfield(r, :pointx)
@inline __getpointy(r::Ray2D{T}) where {T} = getfield(r, :pointy)
@inline __getdirectionx(r::Ray2D{T}) where {T} = getfield(r, :directionx)
@inline __getdirectiony(r::Ray2D{T}) where {T} = getfield(r, :directiony)

@inline __setpointx!(r::Ray2D{T}, value::T) where {T} = setfield!(r, :pointx, value)
@inline __setpointy!(r::Ray2D{T}, value::T) where {T} = setfield!(r, :pointy, value)
@inline __setdirectionx!(r::Ray2D{T}, value::T) where {T} = setfield!(r, :directionx, value)
@inline __setdirectiony!(r::Ray2D{T}, value::T) where {T} = setfield!(r, :directiony, value)
@inline __seti!(r::Ray2D{T}, value::Int) where {T} = setfield!(r, :i, value)
@inline __setj!(r::Ray2D{T}, value::Int) where {T} = setfield!(r, :j, value)

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
  - `pointx::A`: x-coordinate of the intersection point in km.
  - `pointy::A`: y-coordinate of the intersection point in km.
  - `directionx::A`: x-component of the intersection direction vector.
  - `directiony::A`: y-component of the intersection direction vector.
  - `index_i::Ai`: index of the wedge in which the ray intersects.
  - `index_j::Ai`: index of the sub-wedge in which the ray intersects.
  - `length::A`: length of the intersection in km.
  - `altitude::A`: altitude of the intersection in km.
  - `azimuth::A`: azimuth angle of the intersection in degrees.
"""
struct ResultRay{M,T<:AbstractFloat,A<:AbstractVector{T},Ai<:AbstractVector{Int}}
    pointx::A
    pointy::A
    directionx::A
    directiony::A
    index_i::Ai
    index_j::Ai
    length::   A
    altitude:: A
    azimuth::  A

    function ResultRay{T}(n::Int) where {T<:AbstractFloat}
        pointx = zeros(T, n)
        pointy = zeros(T, n)
        directionx = zeros(T, n)
        directiony = zeros(T, n)
        index_i = zeros(Int, n)
        index_j = zeros(Int, n)
        length = zeros(T, n)
        altitude = zeros(T, n)
        azimuth = zeros(T, n)

        return new{n,T,typeof(pointx),typeof(index_i)}(pointx, pointy,
                              directionx, directiony,
                              index_i, index_j,
                              length, altitude,
                              azimuth)
    end

    #ResultRay{T}(directionx::T=0.0, directiony::T=0.0, pointx::T=0.0, pointy::T=0.0,
    #              index_i::Int=0, index_j::Int=0, length::T=0.0,
    #              altitude::T=0.0, azimuth::T=0.0) where T =  new(pointx, pointy, directionx, directiony, index_i, index_j, length, altitude, azimuth)
end

 Base.show(io::IO, ::MIME"text/plain", r::ResultRay{T}) where {T} = begin
    print(io, "ResultRay{", T, "}")
 end

TOEXP=Symbol[]

for prop in (:pointx, :pointy, :directionx, :directiony, :index_i, :index_j, :length, :altitude, :azimuth)
    getfunc = Symbol("get", prop)
    TT      = :T
    docstr_get = """
    `$getfunc(r::ResultRay, int::Int)`
    Return the `:$prop` field of the `ResultRay`.
    """
    @eval begin
      @doc $docstr_get
      @inline function $getfunc(r::ResultRay{M},idx::Int=1) where {M}
        if 1<=idx<=M
        getfield(r, $(QuoteNode(prop)))[idx]
        else
          throw(ArgumentError("Index out of bounds: $idx for ResultRay with length $M"))
        end
      end
    end
    push!(TOEXP, getfunc)
    setfunc = Symbol("__set", prop, "!")
    docstr_set = """
    `$setfunc(r::ResultRay, value)`
    Set the `:$prop` field of the `ResultRay`.
    """
    @eval begin

      @doc $docstr_set
      @inline function $setfunc(r::ResultRay{M,T}, idx::Int, value) where {M,T}
        if 1<=idx<=M
          r.$prop[idx] = value
        else
          throw(ArgumentError("Index out of bounds: $idx for ResultRay with length $M"))
        end
        return r
      end
    end
    push!(TOEXP, setfunc)
end
