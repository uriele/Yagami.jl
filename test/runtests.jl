using Yagami
using Test
using Aqua
using JET
include("MaterialProperties/ciddor_consistency.jl")
include("MaterialProperties/mathar_consistency.jl")

@testset "Yagami.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Yagami)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(Yagami; target_defined_modules = true)
    end
    # Write your tests here.
end
