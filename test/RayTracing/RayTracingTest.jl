using Test
using Yagami.RayTracing
using Yagami.YagamiCore
using GeometryBasics
using Yagami.MaterialProperties: celsius_to_kelvin,atm_to_pascal,kelvin_to_celsius
using Logging,LoggingExtras
using Yagami.RayTracing:__geth, __getθ

@testset "Ray Data Structures" begin
  include("tracingtest.jl")
end
@testset "Atmosphere Interpolation" begin
  include("interpolationtest.jl")
end

@testset "Approximation Accuracy" begin
  include("earthmodeltest.jl")
end

@testset "Create atmosphere" begin
  include("atmospheretest.jl")
end

@info "Testing Cairt File: $testfile"
@info "isfile(testfile)=$(isfile(testfile))"
if isfile(testfile)
  @testset "Cairt File" begin
    include("cairttest.jl")
  end
end
