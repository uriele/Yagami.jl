using Yagami.RayTracing
using Yagami.RayTracing: AbstractResult
using StructArrays
using Test
@testset "Tracing Results Structures" begin

  pointsx= fill(MAJORAXIS(),100)
  pointsy= fill(MAJORAXIS()/cosd(10),100)
  directinsx =zeros(100)
  directinsy = ones(100)
  s1=StructArray{SimpleResult{Float64}}(undef,100,length(pointsx))

  @test eltype(s1) <: AbstractResult


  @testset "SimpleResult Set/Get" begin

    altitude = 120.
    azimuth = 350
    length_t = 1000.
    i = 1
    j = 1
    islevel = true
    desecending = true

    s1.altitude .= altitude
    s1.azimuth .= azimuth
    s1.islevel .= islevel
    s1.descending[:] .= desecending
    s1.i .= i
    s1.j .= j
    s1.length_t .= length_t
    @inbounds for i in eachindex(pointsx)
      s1.pointx[:,i] .= pointsx[i]
      s1.pointy[:,i] .= pointsy[i]
      s1.directionx[:,i] .= directinsx[i]
      s1.directiony[:,i] .= directinsy[i]
    end

    @test all(s1.altitude .== altitude)
    @test all(s1.azimuth .== azimuth)
    @test all(s1.islevel .== islevel)
    @test all(s1.descending .== desecending)
    @test all(s1.i .== i)
    @test all(s1.j .== j)
    @test all(s1.length_t .== length_t)
    for i in eachindex(pointsx)
      @test all(s1.pointx[:,i] .== pointsx[i])
      @test all(s1.pointy[:,i] .== pointsy[i])
      @test all(s1.directionx[:,i] .== directinsx[i])
      @test all(s1.directiony[:,i] .== directinsy[i])
    end
  end
end
