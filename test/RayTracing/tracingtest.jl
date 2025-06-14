using Yagami.RayTracing
using Yagami.RayTracing:__setpointx!, __setpointy!, __setdirectionx!, __setdirectiony!, __seti!, __setj!,__setaltitude!,__setazimuth!,__setlength!
using Yagami.RayTracing:__getaltitude, __getpointx, __getpointy,
__getdirectionx, __getdirectiony, __geti, __getj,__getlength,__getazimuth
using Yagami:MAJORAXIS
using StructArrays

@testset "Tracing Data Structures" begin
  pointsx= fill(MAJORAXIS(),100)
  pointsy= fill(MAJORAXIS()/cosd(10),100)
  directinsx =zeros(100)
  directinsy = ones(100)

  s=StructArray{Ray2D{Float64}, 1}(undef,100)
  __setpointx!(s,1,pointsx[1])


  s1=StructArray{ResultRay{Float64}, 1}(undef,100)
  @simd for i in eachindex(pointsx)
      s.pointx[i] = pointsx[i]
      s.pointy[i] = pointsy[i]
      s.directionx[i] = directinsx[i]
      s.directiony[i] = directinsy[i]
      s.i[i] = 2
      s.j[i] = 1
      __setpointx!(s1,i, pointsx[i])
      __setpointy!(s1,i, pointsy[i])
      __setdirectionx!(s1,i, directinsx[i])
      __setdirectiony!(s1,i, directinsy[i])
      __seti!(s1,i, 2)
      __setj!(s1,i, 1)
      __setaltitude!(s1,i, 10.0)
      __setazimuth!(s1,i, 20.0)
      __setlength!(s1,i, 100.0)
  end

  @test s.pointx== pointsx
  @test s.pointy== pointsy
  @test s.directionx== directinsx
  @test s.directiony== directinsy
  @test s.i == fill(2,100)
  @test s.j == fill(1,100)

  # test helper functions
  @test getpoint(s[1],:x) == getpoint(s[1],1) == s.pointx[1]
  @test getpoint(s[1],:y) == getpoint(s[1],2) == s.pointy[1]
  @test getdirection(s[1],:x) == getdirection(s[1],1) == s.directionx[1]
  @test getdirection(s[1],:y) == getdirection(s[1],2) == s.directiony[1]


  @test getwedgeindex(s[1],:j) == s.j[1] == getwedgeindex(s[1],2)
  @test getwedgeindex(s[1],:i) == s.i[1] == getwedgeindex(s[1],1)
  @test getaltitude(s1[1])==s1.intersection_altitude[1] ==  __getaltitude(s1,1)
  @test getazimuth(s1[1])==s1.intersection_azimuth[1] == __getazimuth(s1,1)
  @test getlength(s1[1])==s1.intersection_length[1] == __getlength(s1,1)
  @test getpoint(s1[1],:x)==s1.intersection_pointx[1] == __getpointx(s1,1)
  @test getpoint(s1[1],:y)==s1.intersection_pointy[1] == __getpointy(s1,1)
  @test getdirection(s1[1],:x)==s1.intersection_directionx[1] == __getdirectionx(s1,1)
  @test getdirection(s1[1],:y)==s1.intersection_directiony[1] == __getdirectiony(s1,1)
  @test getwedgeindex(s1[1],:i)==s1.intersection_index_i[1] == __geti(s1,1)
  @test getwedgeindex(s1[1],:j)==s1.intersection_index_j[1] == __getj(s1,2)
end
