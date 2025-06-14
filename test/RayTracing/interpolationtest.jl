
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
    @test_throws BoundsError __geth(-1.,h,__hmin, __hmax, M)
    @test_throws BoundsError __geth(101.,h,__hmin, __hmax, M)

    @test __getθ(0.01,θ,__θmin, __θmax, N) == N
    @test __getθ(0.0,θ,__θmin, __θmax, N)  == N
    @test __getθ(359.5,θ,__θmin, __θmax, N) == N
    @test __getθ(0.1,θ,__θmin, __θmax, N) == 1

    @test __getθ(θ[1]+0.01,θ,__θmin, __θmax, N) == 1
    @test __getθ(θ[1]-0.01,θ,__θmin, __θmax, N) == N
  end


  # create the interpolation objects

  exact1=[interpfunc1(θ,h) for θ in θ, h in h];
  exact2=[interpfunc2(θ,h) for θ in θ, h in h];
  exact3=[interpfunc3(θ,h) for θ in θ, h in h];
  exact4=[interpfunc4(θ,h) for θ in θ, h in h];

  itp1=AtmInterpolate(θ,h, exact1;logh=true);
  itp2=AtmInterpolate(θ,h, exact2;logh=false);
  itp3=AtmInterpolate(θ,h, exact3;logh=true);
  itp4=AtmInterpolate(θ,h, exact4;logh=true);

  # check the properties of the interpolation objects
  @testset "Test interpolation objects" begin
    @test itp1.interpolation == LogLinear
    @test itp2.interpolation == BiLinear
    @test itp3.interpolation == LogLinear
    @test itp4.interpolation == LogLinear

    @testset "Test bounds" begin
      @test_throws BoundsError itp1(θ[1], h[1]+0.1) # h is out of bounds
      @test_throws BoundsError itp2(θ[1], h[1]+0.1) # h is out of bounds
      @test_throws BoundsError itp3(θ[1], h[1]+0.1) # h is out of bounds
      @test_throws BoundsError itp4(θ[1], h[1]+0.1) # h is out of bounds
      @test_throws BoundsError itp1(θ[1], h[end]-0.1) # h is out of bounds
      @test_throws BoundsError itp2(θ[1], h[end]-0.1) # h is out of bounds
      @test_throws BoundsError itp3(θ[1], h[end]-0.1) # h is out of bounds
      @test_throws BoundsError itp4(θ[1], h[end]-0.1) # h is out of bounds
    end
  end

  @testset "Test Interpolation" begin


    exact1=[interpfunc1(θ,h) for θ in θ, h in h];
    exact2=[interpfunc2(θ,h) for θ in θ, h in h];
    exact3=[interpfunc3(θ,h) for θ in θ, h in h];
    exact4=[interpfunc4(θ,h) for θ in θ, h in h];


    extrap1 = [itp1(θ,h) for θ in θ, h in h];
    extrap2 = [itp2(θ,h) for θ in θ, h in h];
    extrap3 = [itp3(θ,h) for θ in θ, h in h];
    extrap4 = [itp4(θ,h) for θ in θ, h in h];


  # sanity check
    @test exact1== extrap1
    @test exact2== extrap2
    @test exact3== extrap3
    @test exact4== extrap4


    exact1 = [interpfunc1(θ+0.01,h) for θ in θ, h in h];
    exact2 = [interpfunc2(θ+0.01,h) for θ in θ, h in h];
    exact3 = [interpfunc3(θ+0.01,h) for θ in θ, h in h];
    exact4 = [interpfunc4(θ+0.01,h) for θ in θ, h in h];

    extrap1 = [itp1(θ+0.01,h) for θ in θ, h in h];
    extrap2 = [itp2(θ+0.01,h) for θ in θ, h in h];
    extrap3 = [itp3(θ+0.01,h) for θ in θ, h in h];
    extrap4 = [itp4(θ+0.01,h) for θ in θ, h in h];

    # accept a tolerance of 1e-7 for the extrapolated values for a function going from 0 to 100
    @test isapprox(exact1, extrap1,rtol=1e-7)
    @test isapprox(exact2, extrap2,rtol=1e-7)
    @test isapprox(exact3, extrap3,rtol=1e-7)
    @test isapprox(exact4, extrap4,rtol=1e-7)

    exact1 = [interpfunc1(θ+0.01,h-0.1) for θ in θ, h in h[1:end-1]];  # stop 1 shy from the end to stay in bounds
    exact2 = [interpfunc2(θ+0.01,h-0.1) for θ in θ, h in h[1:end-1]];
    exact3 = [interpfunc3(θ+0.01,h-0.1) for θ in θ, h in h[1:end-1]];
    exact4 = [interpfunc4(θ+0.01,h-0.1) for θ in θ, h in h[1:end-1]];

    extrap1 = [itp1(θ+0.01,h-0.1) for θ in θ, h in h[1:end-1]];
    extrap2 = [itp2(θ+0.01,h-0.1) for θ in θ, h in h[1:end-1]];
    extrap3 = [itp3(θ+0.01,h-0.1) for θ in θ, h in h[1:end-1]];
    extrap4 = [itp4(θ+0.01,h-0.1) for θ in θ, h in h[1:end-1]];

    @test isapprox(exact1, extrap1,rtol=1e-7)
    @test isapprox(exact2, extrap2,rtol=1e-7)
    @test isapprox(exact3, extrap3,rtol=1e-6) # this is the sum of the two functions, so it is less accurate
    @test isapprox(exact4, extrap4,rtol=1e-7)
  end
