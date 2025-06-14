module YagamiCore
  using Reexport
  using Unitful: km,m,ustrip,uconvert
  using CoordRefSystems
  using Unitful: km,m, ustrip, uconvert
  @reexport using CoordRefSystems: WGS84Latest
  include("constants.jl")
  include("utilities.jl")


  export MAJORAXIS,MINORAXIS, ECCENTRICITY², ECCENTRICITY, NORMALIZEMINORAXIS, COMPLECCENTRICITY²
  export WGS84MAJORAXIS,WGS84MINORAXIS, WGS84ECCENTRICITY², WGS84ECCENTRICITY, WGS84NORMALIZEMINORAXIS, WGS84COMPLECCENTRICITY²
  export setdatum!,setmajoraxis!,setminoraxis!,seteccentricity²!
end
