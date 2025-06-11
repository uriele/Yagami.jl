const mathar_functions = (
  mathar013_025_refractive_index,
  mathar028_042_refractive_index,
  mathar043_052_refractive_index,
  mathar075_141_refractive_index,
  mathar160_240_refractive_index
)


const mathar_functions! = (
  mathar013_025_refractive_index!,
  mathar028_042_refractive_index!,
  mathar043_052_refractive_index!,
  mathar075_141_refractive_index!,
  mathar160_240_refractive_index!
)

const mathar_pythons=(
  "mathar013_025.py",
  "mathar028_042.py",
  "mathar043_052.py",
  "mathar075_141.py",
  "mathar160_240.py"
)


for (mathar_function,mathar_function!,mathar_python) in zip(mathar_functions,mathar_functions!,mathar_pythons)
  @testset "Consistency of Mathar refractive index with Python implementation for $mathar_function"  begin
    matharpy=load_python_scripts("MaterialProperties/pythonscript/$mathar_python")
    n_py = (λ,t,p,h) -> matharpy.n(pyimport("numpy").array(λ),t,p,h)
    λ= LinRange(0.3,1.69,10) # wavelength in µm
    λ=vcat(λ...,LinRange(9.5,11.0,10)...)
    t= LinRange(-40,100,50) # temperature in °C
    T= @. celsius_to_kelvin(t) # temperature in K
    p= LinRange(800,CIDREF_PRESSUREAIR,50) # pressure in Pa
    h= LinRange(0,1,10) # humidity

    n_yagami = Array{Float64}(undef,length(λ),length(t),length(p),length(h));
    n_python = similar(n_yagami);
    T=collect(T)
    ################################################################################

    for (i,T) in enumerate(T)
      for (j,p) in enumerate(p)
        for (k,h) in enumerate(h)
          @inbounds n_yagami[:,i,j,k]      .= @. mathar_function(T,p,λ,h)
          @inbounds n_python[:,i,j,k]    .= n_py(λ,T,p,h)
        end
      end
    end

    @test n_yagami ≈ n_python

    T! =celsius_to_kelvin.(LinRange(-40,100,100))
    p! = atm_to_pascal.(LinRange(0.001,1,100)) # pressure in Pa
    λ! =10.0*ones(100) # wavelength in µm
    h! =LinRange(0,1,100)


    n_yagami! = similar(T!);
    n_yagami = similar(n_yagami!);
    for (i,(T,p,λ,h)) in enumerate(zip(T!,p!,λ!,h!))
      @inbounds n_yagami[i] = mathar_function(T,p,λ,h)
    end

    mathar_function!(n_yagami!,T!,p!,λ!,h!)

    @test n_yagami ≈ n_yagami!
  end
end
