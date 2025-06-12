using Test
using Yagami.RayTracing
using Yagami.YagamiCore
using GeometryBasics

using Yagami.RayTracing:__geth, __getθ
include("testutils.jl")


#@testset "Atmosphere Interpolation" begin

  h=LinRange(100,0,100)
  θ=LinRange(0.1,359,3600)

  @testset "Test __idxh and __idxθ" begin

      (__hmin, __hmax) = extrema(h)
      (__θmin, __θmax) = extrema(θ)
      (N,M) = (length(θ), length(h))

      idxh=collect(1:M)
      idxθ=collect(1:N)

      _idxh=similar(idxh,M)
      _idxθ=similar(idxθ,N)

      map!(θi-> __getθ(θi,θ,__θmin, __θmax, N), _idxθ, θ)
      map!(hi-> __geth(hi,h,__hmin, __hmax, M), _idxh, h)

      # sanity check
      @test _idxθ == idxθ
      @test _idxh == idxh

      # check if it gets the correct indeces
      _idxh=similar(idxh,M-2)
      _idxθ=similar(idxθ,N)

      θthird = mod.(θ[1:end].+vcat((θ[2:end]-θ[1:end-1])/3,(360-θ[end])/3),360)
      hthird = h[2:end-1]+(h[2:end-1]-h[3:end])/3

      map!(θi-> __getθ(θi,θ,__θmin, __θmax, N), _idxθ, θthird)
      map!(hi-> __geth(hi,h,__hmin, __hmax, M), _idxh, hthird)


    @test _idxθ[1:end-1]== idxθ[1:end-1]
    @test _idxθ[end] == idxθ[end]
    @test _idxh     ==  idxh[1:end-2] # since we have one less h value
    @test_throws ArgumentError __geth(-1.,h,__hmin, __hmax, M)
    @test_throws ArgumentError __geth(101.,h,__hmin, __hmax, M)

    @test __getθ(0.01,θ,__θmin, __θmax, N) == N
    @test __getθ(0.0,θ,__θmin, __θmax, N)  == N
    @test __getθ(359.5,θ,__θmin, __θmax, N) == N
    @test __getθ(0.1,θ,__θmin, __θmax, N) == 1

  end
  ext1=[interpfunc1(θ,h) for θ in θ, h in h];
  ext2=[interpfunc2(θ,h) for θ in θ, h in h];
  ext3=[interpfunc3(θ,h) for θ in θ, h in h];
  ext4=[interpfunc4(θ,h) for θ in θ, h in h];

  itp1=AtmInterpolate(θ,h, ext1;logh=true);
  itp2=AtmInterpolate(θ,h, ext2;logh=false);
  itp3=AtmInterpolate(θ,h, ext3;logh=true);
  itp4=AtmInterpolate(θ,h, ext4;logh=true);

  itp1(10,10)


  itpat_knots1 = [itp1(θ,h) for θ in θ, h in h];
  itpat_knots2 = [itp2(θ,h) for θ in θ, h in h];
  itpat_knots3 = [itp3(θ,h) for θ in θ, h in h];
  itpat_knots4 = [itp4(θ,h) for θ in θ, h in h];


end

@testset "Approximation Accuracy" begin
  h=LinRange(-10,800,800)
  θ=LinRange(1,360,3600)
  MAJORAXIS= WGS84MAJORAXIS
  MINORAXIS= WGS84MINORAXIS
  hθ_fukushima=Point2f[]
  hθ_bowring=Point2f[]
  hθ_exact=Point2f[]
  xy=Point2f[]
  for h in h
    for θ in θ
      (X,Y)=ellipsepoint(h,θ)
      push!(xy,Point2f(X,Y))
      push!(hθ_fukushima,ray2_altitudeangle_fukushima_verbose(X,Y,MAJORAXIS,MINORAXIS))
      push!(hθ_bowring,ray2_altitudeangle_bowring_verbose(X,Y,MAJORAXIS,MINORAXIS))
      push!(hθ_exact,(h,θ))
    end
  end

  @test hθ_fukushima ≈ hθ_exact
  @test hθ_bowring ≈ hθ_exact
end
