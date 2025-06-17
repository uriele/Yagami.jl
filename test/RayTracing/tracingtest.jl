using Yagami
using StructArrays: StructArray
using Test
testfile="$(pwd())/RayTracing/_data/cairt.nc"
prob1=RayTracingProblem(testfile;)
res_serial=StructArray(fill(SimpleResult(),200,50));
res_parallel=StructArray(fill(SimpleResult(),200,50));

_pointx=prob1.pointsx[2]
_pointy=prob1.pointsy[2]
_directionx=prob1.directionsx[2]
_directiony=prob1.directionsy[2]


raytracing!(res_serial,prob1)
raytracing_parallel!(res_parallel,prob1)

tangent_h_base= prob1.tangent_h
tangent_θ_base= mod.(prob1.tangent_θ,360)


tangent_h_serial= similar(tangent_h_base)
tangent_θ_serial= similar(tangent_θ_base)
tangent_h_parallel= similar(tangent_h_base)
tangent_θ_parallel= similar(tangent_θ_base)
@testset "Serial vs Parallel Ray Tracing" begin
  for field in propertynames(res_serial)
    @eval @test res_serial.$field == res_parallel.$field
  end
  @inbounds for j in axes(res_serial,2)
    idx_cut_serial                                  = findfirst(res_serial.length_t[:,j] .== 0)
    idx_cut_parallel                                = findfirst(res_parallel.length_t[:,j] .== 0)
    @test idx_cut_serial == idx_cut_parallel
    tangent_h_serial[j],idx_min_serial              = findmin(res_serial[1:idx_cut_serial-1,j].altitude)
    tangent_θ_serial[j]                             = res_serial[idx_min_serial,j].azimuth
    tangent_h_parallel[j],idx_min_parallel  = findmin(res_parallel[1:idx_cut_parallel-1,j].altitude)
    tangent_θ_parallel[j]                           = res_parallel[idx_min_parallel,j].azimuth
    @test idx_min_serial == idx_min_parallel
    @test tangent_h_serial[j] ≈ tangent_h_parallel[j]
    @test tangent_θ_serial[j] ≈ tangent_θ_parallel[j]
  end
  # difference less than 0.5km
  @test maximum(abs.(tangent_h_serial-tangent_h_base))<0.5
  @test tangent_h_serial== tangent_h_parallel
  # difference less than 1 degree
  @test all(extrema(abs.(tangent_θ_serial-tangent_θ_base)).<1)
  @test tangent_θ_serial== tangent_θ_parallel
end
