using LinearAlgebra,StructArrays,FastGaussQuadrature
using StaticArrays: SVector

"""
     $SIGNATURES

A structure to hold quadrature points and weights for numerical integration. Where `N` is the number of quadrature points and `T` is the type of the points and weights.
"""
struct QuadraturePoints{N,T<:AbstractFloat}
    points::SVector{N,T}
    weights::SVector{N,T}
    function QuadraturePoints(points::A,weights::A) where {T<:AbstractFloat,A<:AbstractVector{T}}
      Np= length(points)
      Nw= length(weights)
      if Np != Nw
        throw(ArgumentError("Points and weights must have the same length."))
      end
      N=Np
      new{N,T}(SVector{N,T}(points), SVector{N,T}(weights))
    end
end


"""
  $SIGNATURES
Get quadrature points and weights for numerical integration using Gauss-Legendre quadrature.
  The function returns a `QuadraturePoints` object containing the points and weights.

  # Arguments
  - `T`: Type of the quadrature points (default is `Float64`).
  - `n`: Number of quadrature points (must be positive).

  # Returns
  A `QuadraturePoints` object with the quadrature points and weights.
"""
function get_quadrature_points(::Type{T},n::Int) where {T<:AbstractFloat}
    if n <= 0
        throw(ArgumentError("Number of quadrature points must be positive."))
    end
    points, weights = gausslegendre(n)
    # my domain will alwaus be between 0 and t
    return QuadraturePoints(points, weights)
end
get_quadrature_points(n::Int) = get_quadrature_points(Float64,n)



"""
     $SIGNATURES
Compute the linear integral of a function `f` over the quadrature points `Q` scaled by a factor `t`. It is assumed that the integral starts at `0` and ends at `t`.
The function uses the quadrature points and weights to evaluate the integral numerically.
The function evaluates `f` at the quadrature points scaled by `t`, computes the dot product with the weights, and returns the result multiplied by `t/2`.

The function f needs to have the following signature:
```julia
f(points::AbstractVector{T}, p) where {T<:AbstractFloat}
```
where p is an optional parameter that can be passed to the function to account for additional parameters in the integral.

  # Arguments
  - `f`: A function that takes a vector of points and an optional parameter `p`.
  - `Q`: A `QuadraturePoints` object containing the quadrature points and weights.
  - `t`: the end of the integral path, which is a scalar value.
  - `p`: An optional parameter passed to the function `f`.

  # Returns
  The computed linear integral as a scalar value.
"""
function linintegral(f::F,Q::QuadraturePoints{N,T},t::T=T(1.0),p=nothing) where {F,T<:AbstractFloat,N}
    t_half=t*convert(T,0.5)
    dot(Q.weights, f((@. Q.points*t_half+t_half),p)) * t_half
end



"""
     $SIGNATURES
Compute the linear integral of a function `f` over the quadrature points `Q` scaled by a factor `t`. It is assumed that the integral starts at `0` and ends at `t`.
The function uses the quadrature points and weights to evaluate the integral numerically.
The function evaluates `f` at the quadrature points scaled by `t`, computes the dot product with the weights, and stores the result in the output array `out`.
The function f needs to have the following signature:
```julia
f(points::AbstractVector{T}, p) where {T<:AbstractFloat}
```
where p is an optional parameter that can be passed to the function to account for additional parameters in the integral.

  # Arguments
  - `out`: An output array to store the results of the integral.
  - `f`: A function that takes a vector of points and an optional parameter `p`.
  - `Q`: A `QuadraturePoints` object containing the quadrature points and weights.
  - `t`: the end of the integral path, which is a vector of the same size as `out`.
  - `p`: An optional parameter passed to the function `f`. Can be a vector of the same size of `out` or a single value.

  # Returns
  The output array with the computed linear integrals.
"""
function linintegral!(out::AbstractArray{T},f::F,Q::QuadraturePoints{N,T},t::AbstractArray{T},p::AbstractArray) where {T<:AbstractFloat,N,F}
    @simd  for i in eachindex(out)
        @inbounds t_half=t[i]*convert(T,0.5)
        @inbounds out[i]=dot(Q.weights, f((@.Q.points*t_half+t_half),p[i])) * t_half
    end
    return out
end

function linintegral!(out::AbstractArray{T},f::F,Q::QuadraturePoints{N,T},t::AbstractArray{T},p::S) where {T<:AbstractFloat,N,S,F}
    @simd  for i in eachindex(out)
        @inbounds t_half=t[i]*convert(T,0.5)
        @inbounds out[i]=dot(Q.weights, f((@.Q.points*t_half+t_half),p)) * t_half
    end
    return out
end
