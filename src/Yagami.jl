module Yagami
using Reexport
include("YagamiCore/YagamiCore.jl")
include("MaterialProperties/MaterialProperties.jl")
include("RayTracing/RayTracing.jl")
include("CurtisGodson/CurtisGodson.jl")
# Write your package code here.
#using .YagamiCore
#using .MaterialProperties
#using .RayTracing
#using .CurtisGodson


@reexport using .YagamiCore: MAJORAXIS,MINORAXIS, ECCENTRICITY², ECCENTRICITY, NORMALIZEMINORAXIS, COMPLECCENTRICITY²
@reexport using .YagamiCore: WGS84MAJORAXIS,WGS84MINORAXIS, WGS84ECCENTRICITY², WGS84ECCENTRICITY, WGS84NORMALIZEMINORAXIS, WGS84COMPLECCENTRICITY²
@reexport using .YagamiCore: setdatum!,setmajoraxis!,setminoraxis!,seteccentricity²!


@reexport using .MaterialProperties: refractive_index, refractive_index!
@reexport using .RayTracing: snellslaw!

@reexport using .MaterialProperties: AirModel, Ciddor, Mathar, Carlotti
@reexport using .MaterialProperties: Mathar160_240, Mathar075_141, Mathar043_052, Mathar028_042, Mathar013_025
@reexport using .RayTracing: gettemperature, getpressure, gethumidity, getco2ppm

@reexport using .RayTracing: AtmInterpolate, BiLinear, LogLinear
@reexport using .RayTracing: MeanType, ArithmeticMean, GeometricMean, LogMean

@reexport using .RayTracing: create_atmosphere, create_hlevelset,create_radii

@reexport using .RayTracing: Ray2D, ResultRay
@reexport using .RayTracing: Ray2Ds, ResultRays
@reexport using .RayTracing: getwedgeindex, getpoint, getdirection,getaltitude,getazimuth,getlength
@reexport using .RayTracing: ray2_angle, ray2_altitude, ray2_altitudeangle
@reexport using .RayTracing: EarthApproximation, Fukushima, Bowring
@reexport using .RayTracing: ray2_altitudeangle_bowring_verbose, ray2_altitudeangle_fukushima_verbose,
ray2_angle_bowring_verbose, ray2_angle_fukushima_verbose,
ray2_altitude_bowring_verbose, ray2_altitude_fukushima_verbose
@reexport using .RayTracing: ray2_altitudeangle_bowring, ray2_altitudeangle_fukushima,
ray2_angle_bowring, ray2_angle_fukushima,
ray2_altitude_bowring, ray2_altitude_fukushima

end
