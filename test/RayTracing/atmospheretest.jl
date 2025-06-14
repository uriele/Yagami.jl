using Test
using Yagami
using Yagami.MaterialProperties: celsius_to_kelvin, atm_to_pascal, kelvin_to_celsius,pascal_to_atm
using Yagami.RayTracing: AtmosphereSetting, create_atmosphere
using Yagami.RayTracing:grid_refractiveindex
using Yagami.RayTracing

h=create_hlevelset()
θ=create_radii()

@testset "Atmosphere Interpolation" begin
  temperature = [temp_func(θ, h) for θ in θ, h in h]
  pressure = [press_func(θ, h) for θ in θ, h in h]
  co2ppm = fill(0.0, length(θ), length(h))
  humidity = fill(0.0, length(θ), length(h))


  atm1=AtmosphereSetting(θ,h, temperature,pressure,0.0,0.0)
  atm2=AtmosphereSetting(θ,h, temperature,pressure,humidity,0.0)
  atm3=AtmosphereSetting(θ,h, temperature,pressure,humidity,co2ppm)
  @testset "test atmosphere creation" begin
    @test gettemperature(atm1, 10.0,10.0)==gettemperature(atm2, 10.0,10.0)==gettemperature(atm3, 10.0,10.0)
    @test getpressure(atm1, 10.0,10.0)==getpressure(atm2, 10.0,10.0)==getpressure(atm3, 10.0,10.0)
    @test gethumidity(atm1, 10.0,10.0)==gethumidity(atm2, 10.0,10.0)==gethumidity(atm3, 10.0,10.0)
    @test getco2ppm(atm1, 10.0,10.0)==getco2ppm(atm2, 10.0,10.0)==getco2ppm(atm3, 10.0,10.0)
  end


    # Syntetic atmosphere functions for testing
  θᵢ = create_radii(Float64;)
  hᵢ = create_hlevelset(Float64;hmin=4.0)
  temperatureᵢ = Matrix([temp_func(θ, h) for θ in θᵢ, h in hᵢ])
  pressureᵢ = Matrix([press_func(θ, h) for θ in θᵢ, h in hᵢ])
  knots_θ = create_radii(Float64;θmin=0.0, θmax=350.0, radii=36)
  knots_h = create_hlevelset(Float64;hmin=4.0, hmax=120.0, levels=20)



    atm=create_atmosphere(;θᵢ=θᵢ,
      hᵢ=hᵢ,
      temperatureᵢ=temperatureᵢ,
      pressureᵢ=pressureᵢ,
      knots_θ=knots_θ,
      knots_h=knots_h,
      humidity=0.0,
      co2ppm=0.0,
      wavelength=10.0)

    _temperatureᵢ= Matrix{Float64}(undef, length(θᵢ), length(hᵢ))
    _pressureᵢ = Matrix{Float64}(undef, length(θᵢ), length(hᵢ))

    @simd for j in eachindex(hᵢ)
      @inbounds for i in eachindex(θᵢ)
        _temperatureᵢ[i,j]=gettemperature(atm, θᵢ[i], hᵢ[j])
        _pressureᵢ[i,j]=getpressure(atm, θᵢ[i], hᵢ[j])
      end
    end

    @inline mse(x, y) = sum((x .- y).^2) / length(x)

    @testset "test temperature and pressure average" begin

      # sind is a complex function so I allow some error in downscaling
      @test mse(_temperatureᵢ, temperatureᵢ)<5
      @test isapprox(_pressureᵢ, pressureᵢ)
    end
    n =grid_refractiveindex(atm;) # just to check if it works
    n1=grid_refractiveindex(atm; model=Ciddor())
    n2=grid_refractiveindex(atm; meantype=GeometricMean())
    @testset "create atmosphere" begin

      @test n==n1==n2
      n_mathar=grid_refractiveindex(atm; model=Mathar075_141())
      n_ciddorlog=grid_refractiveindex(atm; meantype=LogMean())

      @test n_mathar≠n_ciddorlog # Mathar and Ciddor log mean are not the same
  end
end

@testset "Snell's Law" begin
  normal= [-1.0,0.0] # consider a outward normal to the surface
  direction= [sind(45.0), cosd(45.0)]
  direction_test  = copy(direction)

  snellslaw!(;normal=normal, direction=direction)
  @test direction==direction_test # direction should not change

  # sin(θ1)/sin(θ2) = n2/n1
  sindout=cosd(45.)/1.5
  direction_test2 = [sqrt(1-sindout^2), sindout]

  snellslaw!(;normal=normal, direction=direction,n_transmitted=1.5)
  @test direction ≈ direction_test2
end
