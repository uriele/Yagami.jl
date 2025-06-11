module RayTracing
using LinearAlgebra, StaticArrays
using DocStringExtensions
include("../YagamiCore/YagamiCore.jl")
include("../MaterialProperties/MaterialProperties.jl")
using .YagamiCore
using .MaterialProperties
include("grid/atmospheresetting.jl")
include("constants.jl")
include("models/fukushima.jl")
include("models/bowring.jl")
include("tracing/fastmarchingintersection.jl")

export create_atmosphere, create_hlevelset,create_radii
for model in EXISTINGMODELS
  for what in RETURNWHAT
    funcname=Symbol("ray2",what,"_",model)
    funcname_verbose=Symbol(funcname,"_verbose")
    @eval export $funcname, $funcname_verbose
  end
end


end
