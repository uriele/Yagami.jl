module RayTracing
using LinearAlgebra, StaticArrays
using DocStringExtensions
include("../YagamiCore/YagamiCore.jl")
include("../MaterialProperties/MaterialProperties.jl")
using .YagamiCore
using .MaterialProperties
import NaturalNeighbours.interpolate as NNinterpolate

include("constants.jl")
include("models/fukushima.jl")
include("models/bowring.jl")
include("grid/interpolation.jl")
include("grid/atmospheresetting.jl")
include("tracing/snellslaw.jl")
include("tracing/fastmarchingintersection.jl")
include("utilities.jl")


export snellslaw
export AirModel, Ciddor, Mathar, Carlotti
export Mathar160_240, Mathar075_141, Mathar043_052, Mathar028_042, Mathar013_025
export gettemperature, getpressure, gethumidity, getco2ppm

export AtmInterpolate, BiLinear, LogLinear
export MeanType, ArithmeticMean, GeometricMean, LogMean

export create_atmosphere, create_hlevelset,create_radii
for model in EXISTINGMODELS
  for what in RETURNWHAT
    funcname=Symbol("ray2",what,"_",model)
    funcname_verbose=Symbol(funcname,"_verbose")
    @eval export $funcname, $funcname_verbose
  end
end

export create_hlevelset, create_radii,create_atmosphere,grid_refractiveindex

end
