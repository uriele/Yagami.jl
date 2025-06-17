using Yagami
using Test
using Aqua
using JET
include("testutils.jl")

@testset "Yagami/MaterialProperties" begin
    include("MaterialProperties/MaterialPropertiesTest.jl")
end

@testset "Yagami/RayTracing" begin
  include("RayTracing/RayTracingTest.jl")
end

@testset "Yagami.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Yagami; stale_deps=false)
    end
    @testset "Code linting (JET.jl)" begin
      #  JET.test_package(Yagami; target_defined_modules = true)
    end
    # Write your tests here.
end
