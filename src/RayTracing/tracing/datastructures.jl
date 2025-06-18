struct NoFunc end
abstract type AbstractResult{T<:AbstractFloat} end



struct SimpleResult{T} <:AbstractResult{T}
  pointx::T
  pointy::T
  directionx::T
  directiony::T
  altitude::T
  azimuth::T
  length_t::T
  i::Int
  j::Int
  descending::Bool
  islevel::Bool
end


"""
   $SIGNATURES

Create a `SimpleResult` with default values for the fields. SimpleResult is a
data structure used to store the results of ray tracing calculations and it is a sybtype of
`AbstractResult`.

  # Arguments
- `pointx::T=MAJORAXIS()`: X-coordinate of the point, default is the major axis.
- `pointy::T=MINORAXIS()`: Y-coordinate of the point, default is the minor axis.
- `directionx::T=0.0`: X-component of the direction vector, default is 0.0.
- `directiony::T=-1.0`: Y-component of the direction vector, default is -1.0.
- `altitude::T=zero(T)`: Altitude, default is 0.
- `azimuth::T=zero(T)`: Azimuth, default is 0.
- `length::T=zero(T)`: Length, default is 0.
- `i::Int=0`: Index i, default is 0.
- `j::Int=0`: Index j, default is 0.
- `descending::Bool=true`: Whether the ray is descending, default is true.

  # Returns
- A `SimpleResult` object with the specified or default values.

"""
SimpleResult(pointx::T=MAJORAXIS(), pointy::T=MINORAXIS(), directionx::T=0.0, directiony::T=-1.0,
  altitude::T=zero(T), azimuth::T=zero(T), length::T=zero(T), i::Int=0, j::Int=0,descending::Bool=true) where {T<:AbstractFloat} =
  SimpleResult{T}(pointx, pointy, directionx, directiony, altitude, azimuth, length, i, j,descending,true)


  azimuthlocal(results::AR,i) where {AR<:AbstractMatrix{<:AbstractResult}} = begin
    idx= findfirst(results.length_t[:,i] .== 0)
    results.azimuth[1:idx-1,i]
  end
  altitudelocal(results::AR,i) where {AR<:AbstractMatrix{<:AbstractResult}} = begin
    idx= findfirst(results.length_t[:,i] .== 0)
    results.altitude[1:idx-1,i]
  end
  pointxlocal(results::AR,i) where {AR<:AbstractMatrix{<:AbstractResult}} = begin
    idx= findfirst(results.length_t[:,i] .== 0)
    results.pointx[1:idx-1,i]
  end
  pointylocal(results::AR,i) where {AR<:AbstractMatrix{<:AbstractResult}} = begin
    idx= findfirst(results.length_t[:,i] .== 0)
    results.pointy[1:idx-1,i]
  end

  lengthlocal(results::AR,i) where {AR<:AbstractMatrix{<:AbstractResult}} = begin
    idx= findfirst(results.length_t[:,i] .== 0)
    results.length_t[1:idx-1,i]
  end
