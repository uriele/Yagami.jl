using Test
using Yagami.MaterialProperties: CIDREF_PRESSUREAIR
using Yagami.MaterialProperties:celsius_to_kelvin, atm_to_pascal
using Yagami.MaterialProperties:__saturation_vapor_pressure
using Yagami
using PyCall
using Yagami.MaterialProperties


include("ciddor_consistency.jl")
include("mathar_consistency.jl")

# Test interfacest

t= LinRange(-40,100,100) # temperature in °C
T= @. celsius_to_kelvin(t) # temperature in K
patm= LinRange(0.001,1,100) # pressure in atm
p   = @. atm_to_pascal(patm) # pressure in Pa
h   = LinRange(0,1,100)
co2ppm= LinRange(0,500,100) # CO2 concentration in ppm
λ     = 10.0*ones(100) # wavelength in µm
n     = similar(λ)

@testset "Test intefraces" begin
  @testset "Ciddor interfaces" begin
    @test refractive_index(Ciddor(), T[1], p[1], λ[1], h[1], co2ppm[1]) == ciddor_refractive_index(T[1], p[1], λ[1], h[1], co2ppm[1])
    @test refractive_index!(Ciddor(), n, T, p, λ, h, co2ppm)==  ciddor_refractive_index!(n, T, p, λ, h,co2ppm)
  end
  @testset "Mathar interfaces" begin
    models = (
      Mathar160_240(),
      Mathar075_141(),
      Mathar043_052(),
      Mathar028_042(),
      Mathar013_025()
    )
    modelfuncs = (
      mathar160_240_refractive_index,
      mathar075_141_refractive_index,
      mathar043_052_refractive_index,
      mathar028_042_refractive_index,
      mathar013_025_refractive_index
    )
    modelfuncs! = (
      mathar160_240_refractive_index!,
      mathar075_141_refractive_index!,
      mathar043_052_refractive_index!,
      mathar028_042_refractive_index!,
      mathar013_025_refractive_index!
    )


    for i in eachindex(models)
      @testset "Mathar $(models[i]) interfaces" begin

       model=models[i]
        modelfunc=modelfuncs[i]
        modelfunc! = modelfuncs![i]

        @test refractive_index(model,T[1],p[1],λ[1],h[1],co2ppm[1])==  modelfunc(T[1], p[1], λ[1], h[1])
        @test refractive_index!(model, n, T, p, λ, h, co2ppm)==modelfunc!(n, T, p, λ, h)


        _h      = fill(0.5,100)
        _co2ppm = fill(450.0,100) # CO2 concentration in ppm
        _λ      = fill(10.0,100) # wavelength in µm


        @test refractive_index!(model,n,T,p,_λ[1],_h[1],_co2ppm[1])== modelfunc!(n, T, p, λ[1], h[1])
        @test refractive_index!(model,n,T,p,_λ,_h,_co2ppm)== refractive_index!(model,n,T,p,_λ[1],_h[1])
      end
    end
  end
end
