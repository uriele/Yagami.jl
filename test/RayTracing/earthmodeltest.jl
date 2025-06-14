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
