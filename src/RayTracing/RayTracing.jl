module RayTracing
using LinearAlgebra, StaticArrays
using DocStringExtensions
using StructArrays
using Reexport
include("../YagamiCore/YagamiCore.jl")
include("../MaterialProperties/MaterialProperties.jl")
using .YagamiCore
using .MaterialProperties
@reexport using ..MaterialProperties: Mathar160_240, Mathar075_141, Mathar043_052, Mathar028_042, Mathar013_025
@reexport using ..MaterialProperties: Ciddor, Mathar, Carlotti, AirModel
@reexport using ..MaterialProperties: refractive_index, refractive_index!
import NaturalNeighbours.interpolate as NNinterpolate

include("constants.jl")

abstract type EarthApproximation end
struct Fukushima <: EarthApproximation end
struct Bowring <: EarthApproximation end


include("models/fukushima.jl")
include("models/bowring.jl")
include("grid/interpolation.jl")
include("grid/atmospheresetting.jl")
include("tracing/datastructures.jl")
include("tracing/snellslaw.jl")
include("tracing/fastmarchingintersection.jl")
include("utilities.jl")

export EarthApproximation, Fukushima, Bowring
ray2_angle(::EA,W::T,Z::T,a::T,b::T) where {T<:AbstractFloat,EA<:EarthApproximation} = throw(ArgumentError("Unknown Earth Approximation: $EA. Use Fukushima or Bowring."))
ray2_altitude(::EA,W::T,Z::T,a::T,b::T) where {T<:AbstractFloat,EA<:EarthApproximation} = throw(ArgumentError("Unknown Earth Approximation: $EA. Use Fukushima or Bowring."))
ray2_altitudeangle(::EA,W::T,Z::T,a::T,b::T) where {T<:AbstractFloat,EA<:EarthApproximation} = throw(ArgumentError("Unknown Earth Approximation: $EA. Use Fukushima or Bowring."))


export Ray2D,ResultRay
export Ray2Ds,ResultRays
export getwedgeindex, getpoint, getdirection,getaltitude,getazimuth,getlength


export snellslaw!
export gettemperature, getpressure, gethumidity, getco2ppm

export AtmInterpolate, BiLinear, LogLinear
export MeanType, ArithmeticMean, GeometricMean, LogMean

export create_atmosphere, create_hlevelset,create_radii
for model in EXISTINGMODELS
  for what in RETURNWHAT
    interfunname=Symbol("ray2",what)
    funcname=Symbol("ray2",what,"_",model)
    funcname_verbose=Symbol(funcname,"_verbose")
    @info "Exporting function $funcname, $funcname_verbose, $interfunname"
    @eval export $funcname, $funcname_verbose, $interfunname
  end
end

export create_hlevelset, create_radii,create_atmosphere,grid_refractiveindex


end
