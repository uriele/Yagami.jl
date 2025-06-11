
include("testutils.jl")

@testset "Consistency of Ciddor refractive index with Python implementation" begin

  ciddorpy=load_python_scripts("MaterialProperties/pythonscript/ciddor.py")
  n_py = (λ,t,p,h,co2ppm) -> ciddorpy.n(pyimport("numpy").array(λ),t,p,h,co2ppm)
  svp_py = (t)-> ciddorpy.saturation_vapor_pressure(pyimport("numpy").array(t))
  Z_py = (t,p,xw) -> ciddorpy.Z(pyimport("numpy").array(t,p,xw))
  λ= LinRange(0.3,1.69,10) # wavelength in µm
  λ=vcat(λ...,LinRange(9.5,11.0,10)...)
  t= LinRange(-40,100,50) # temperature in °C
  T= @. celsius_to_kelvin(t) # temperature in K
  p= LinRange(800,CIDREF_PRESSUREAIR,50) # pressure in Pa
  h= LinRange(0,1,10) # humidity in fraction
  co2ppm= LinRange(0,500,20) # CO2 concentration in ppm

  n_yagami = Array{Float64}(undef,length(λ),length(t),length(p),length(h),length(co2ppm));
  n_python = similar(n_yagami);
  svp_yagami = similar(n_yagami);
  svp_python = similar(n_yagami);

  # correct for really small values since the equation is not continuous across 0
  T=collect(T)
  ################################################################################
  @test(__saturation_vapor_pressure.(T)≈svp_py(T))

  for (i,(t,T)) in enumerate(zip(t,T))
    for (j,p) in enumerate(p)
      for (k,h) in enumerate(h)
        for (l,co2ppm) in enumerate(co2ppm)
            @inbounds n_yagami[:,i,j,k,l]      .= @. Yagami.MaterialProperties.ciddor_refractive_index(T,p,λ,h,co2ppm)
            @inbounds n_python[:,i,j,k,l]    .= n_py(λ,T,p,h,co2ppm)
            @inbounds svp_yagami[:,i,j,k,l]    .=@. Yagami.MaterialProperties.__saturation_vapor_pressure(T)
            @inbounds svp_python[:,i,j,k,l]  .= svp_py(T)
        end
      end
    end
  end

  n_py(λ[1],T[1],p[1],h[1],400.)[] # test single value
  ciddor_refractive_index(T[1],p[1],λ[1],h[1],400.) # test single value

  @test svp_yagami ≈ svp_python
  @test n_yagami ≈ n_python


  T! =celsius_to_kelvin.(LinRange(-40,100,100))
  p! = atm_to_pascal.(LinRange(0.001,1,100)) # pressure in Pa
  λ! =10.0*ones(100) # wavelength in µm
  h! =LinRange(0,1,100)
  co2ppm! = 400.0*ones(25) # CO2 concentration in ppm
  co2ppm! = vcat(co2ppm!,425.0*ones(50)) # CO2 concentration in ppm
  co2ppm! = vcat(co2ppm!,450.0*ones(25)) # CO2 concentration in ppm

  n_yagami! = similar(T!);
  n_yagami = similar(n_yagami!);
  for (i,(T,p,λ,h,co2ppm)) in enumerate(zip(T!,p!,λ!,h!,co2ppm!))
    @inbounds n_yagami[i] = ciddor_refractive_index(T,p,λ,h,co2ppm)
  end

  ciddor_refractive_index!(n_yagami!,T!,p!,λ!,h!,co2ppm!)

  @test n_yagami ≈ n_yagami!

end
