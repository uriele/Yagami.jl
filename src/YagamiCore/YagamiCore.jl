module YagamiCore
  using Reexport
  using Logging, LoggingExtras
  using Unitful: km,m,ustrip,uconvert
  using CoordRefSystems
  using Printf: @sprintf
  using Unitful: km,m, ustrip, uconvert
  @reexport using CoordRefSystems: WGS84Latest
  include("constants.jl")
  include("utilities.jl")
  include("zbrent.jl")

  export MAJORAXIS,MINORAXIS, ECCENTRICITY², ECCENTRICITY, NORMALIZEMINORAXIS, COMPLECCENTRICITY²
  export WGS84MAJORAXIS,WGS84MINORAXIS, WGS84ECCENTRICITY², WGS84ECCENTRICITY, WGS84NORMALIZEMINORAXIS, WGS84COMPLECCENTRICITY²
  export setdatum!,setmajoraxis!,setminoraxis!,seteccentricity²!

  export Zbrent, findraymin, bracketmin
  export clampangle,intersectionrayray

  # find angle from nadir and limb convention
  export nadir_angle_normal, limb_angle_normal
  export infologger
  export file_to_array
end
