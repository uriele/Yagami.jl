using Yagami
using BenchmarkTools

SUITE = BenchmarkGroup()
SUITE["rand"] = @benchmarkable rand(10)

# Write your benchmarks here.
SUITE["create_hlevelset"] = @benchmarkable create_hlevelset(Float64; hmin=4.0, hmax=120.0, levels=20)
SUITE["create_radii"] = @benchmarkable create_radii(Float64; θmin=0.0, θmax=350.0, radii=36)


SUITE["create_atmosphere"] = @benchmarkable create_atmosphere(
    θᵢ=create_radii(Float64; θmin=0.0, θmax=350.0, radii=36),
    hᵢ=create_hlevelset(Float64; hmin=4.0, hmax=120.0, levels=10),
    temperatureᵢ=Matrix{Float64}(undef, 36, 20),
    pressureᵢ=Matrix{Float64}(undef, 36, 20),
    knots_θ=create_radii(Float64; θmin=0.0, θmax=350.0, radii=360),
    knots_h=create_hlevelset(Float64; hmin=4.0, hmax=120.0, levels=20),
    humidity=0.0,
    co2ppm=0.0,
    wavelength=10.0
)
