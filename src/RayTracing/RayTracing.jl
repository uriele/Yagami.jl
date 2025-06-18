module RayTracing
using LinearAlgebra

using StaticArrays:SVector
using Logging,LoggingExtras
using DocStringExtensions
using StructArrays
using Lazy: @forward  # to forward methods of DistanceFunc to Zbrent
using Moshi.Match: @match
using Polyester: @batch
using Reexport
include("../YagamiCore/YagamiCore.jl")
include("../MaterialProperties/MaterialProperties.jl")
using .YagamiCore
using .YagamiCore: infologger
using .MaterialProperties
@reexport using ..MaterialProperties: Mathar160_240, Mathar075_141, Mathar043_052, Mathar028_042, Mathar013_025
@reexport using ..MaterialProperties: Ciddor, Mathar, Carlotti, AirModel
@reexport using ..MaterialProperties: refractive_index, refractive_index!

@reexport using ..YagamiCore: bracketmin, findraymin, Zbrent,__setbracket!
import NaturalNeighbours.interpolate as NNinterpolate

include("constants.jl")

abstract type EarthApproximation end
struct Fukushima <: EarthApproximation end
struct Bowring <: EarthApproximation end

include("models/fukushima.jl")
include("models/bowring.jl")
include("grid/interpolation.jl")
include("grid/atmospheresetting.jl")
include("utilities.jl")
include("tracing/datastructures.jl")
include("tracing/distancefunc.jl")
include("files/cairt.jl")
include("files/loggersimulation.jl")
include("tracing/snellslaw.jl")
include("tracing/minimizationintersection.jl")


export RayTracingProblem,AtmosphereSetting
export EarthApproximation, Fukushima, Bowring
export earthmodelfunction
export DistanceFunc
ray2_angle(::EA,W::T,Z::T,a::T,b::T) where {T<:AbstractFloat,EA<:EarthApproximation} = throw(ArgumentError("Unknown Earth Approximation: $EA. Use Fukushima or Bowring."))
ray2_altitude(::EA,W::T,Z::T,a::T,b::T) where {T<:AbstractFloat,EA<:EarthApproximation} = throw(ArgumentError("Unknown Earth Approximation: $EA. Use Fukushima or Bowring."))
ray2_altitudeangle(::EA,W::T,Z::T,a::T,b::T) where {T<:AbstractFloat,EA<:EarthApproximation} = throw(ArgumentError("Unknown Earth Approximation: $EA. Use Fukushima or Bowring."))


export  getpoint, getdirection


export snellslaw!
export gettemperature, getpressure, gethumidity, getco2ppm
export getknotsh, getknotsθ

export AtmInterpolate, BiLinear, LogLinear
export MeanType, ArithmeticMean, GeometricMean, LogMean

export create_atmosphere, create_hlevelset,create_radii, grid_refractiveindex
for model in EXISTINGMODELS
  for what in RETURNWHAT
    interfunname=Symbol("ray2",what)
    funcname=Symbol("ray2",what,"_",model)
    funcname_verbose=Symbol(funcname,"_verbose")

    @eval export $funcname, $funcname_verbose, $interfunname
  end
end


for func in EXPCAIRTHELPERS
  @eval export $func
end
export create_hlevelset, create_radii,create_atmosphere,grid_refractiveindex


# Main functions for ray tracing
export raytracing!, raytracing_parallel!

export SimpleResult,AbstractResult

export azimuthlocal, altitudelocal, pointxlocal, pointylocal, lengthlocal

end
