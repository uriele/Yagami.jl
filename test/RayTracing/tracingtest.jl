using Yagami.RayTracing
using Yagami.RayTracing:__setpointx!, __setpointy!, __setdirectionx!, __setdirectiony!, __seti!, __setj!,__setaltitude!,__setazimuth!,__setlength!
using Yagami.RayTracing: __getpointx, __getpointy,__getdirectionx, __getdirectiony
using Yagami:MAJORAXIS
using StructArrays
using Test
@testset "Tracing Data Structures" begin
  pointsx= fill(MAJORAXIS(),100)
  pointsy= fill(MAJORAXIS()/cosd(10),100)
  directinsx =zeros(100)
  directinsy = ones(100)

  s = StructArray{Ray2D{Float64}}(undef,100)
  for i in eachindex(s)
    s.pointx[i] = pointsx[i]
    s.pointy[i] = pointsy[i]
    s.directionx[i] = directinsx[i]
    s.directiony[i] = directinsy[i]
    s.i[i] = 1
    s.j[i] = 2
  end

  rex=ResultRay{Float64}(10)
  s1=StructArray(fill(rex,100))
  for i in 1:10
    __setpointx!(s1[1],i,round(i))
    __setpointy!(s1[1],i,round(i+1))
    __setdirectionx!(s1[1],i,round(i+2))
    __setdirectiony!(s1[1],i,round(i+3))
    __setindex_i!(s1[1],i,i)
    __setindex_j!(s1[1],i,i)
    __setaltitude!(s1[1],i,round(i+6))
    __setazimuth!(s1[1],i,round(i+7))
    __setlength!(s1[1],i,round(i+8))
  end

  # test helper functions
  @test getpoint(s[1],:x) == getpoint(s[1],1) == s.pointx[1]
  @test getpoint(s[1],:y) == getpoint(s[1],2) == s.pointy[1]
  @test getdirection(s[1],:x) == getdirection(s[1],1) == s.directionx[1]
  @test getdirection(s[1],:y) == getdirection(s[1],2) == s.directiony[1]

  @test getindex_i(s1[1]) == s1[1].index_i[1] == getindex_j(s1[1],1)
  @test s1[1].index_j[4] == getindex_j(s1[1],4)
  @test s1[1].length[1] == getlength(s1[1],1)
  @test s1[1].altitude[1] == getaltitude(s1[1],1)
  @info s1[1].azimuth[1] == getazimuth(s1[1],1)
  @test getpointx(s1[1])== s1[1].pointx[1]
  @test getpointy(s1[1]) == s1[1].pointy[1]
  @test getdirectionx(s1[1],2) ==  s1[1].directionx[2]
  @test getdirectiony(s1[1],3) ==  s1[1].directiony[3]
end
