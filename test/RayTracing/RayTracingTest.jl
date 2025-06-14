using Test
using Yagami.RayTracing
using Yagami.YagamiCore
using GeometryBasics
using Yagami.MaterialProperties: celsius_to_kelvin,atm_to_pascal,kelvin_to_celsius

using Yagami.RayTracing:__geth, __getθ
include("testutils.jl")


@testset "Atmosphere Interpolation" begin
  include("interpolationtest.jl")
end

@testset "Approximation Accuracy" begin
  include("earthmodeltest.jl")
end

@testset "Create atmosphere" begin
  include("atmospheretest.jl")
end
