module Yagami
include("YagamiCore/YagamiCore.jl")
include("MaterialProperties/MaterialProperties.jl")
include("RayTracing/RayTracing.jl")
include("CurtisGodson/CurtisGodson.jl")
# Write your package code here.
using .YagamiCore

  export MAJORAXIS,MINORAXIS, ECCENTRICITY², ECCENTRICITY, NORMALIZEMINORAXIS, COMPLECCENTRICITY²
  export WGS84MAJORAXIS,WGS84MINORAXIS, WGS84ECCENTRICITY², WGS84ECCENTRICITY, WGS84NORMALIZEMINORAXIS, WGS84COMPLECCENTRICITY²
  export setdatum!,setmajoraxis!,setminoraxis!,seteccentricity²!,seteccentricity!,setnormalizeminoraxis!
end
