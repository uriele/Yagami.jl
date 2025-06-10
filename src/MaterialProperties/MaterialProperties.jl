
module MaterialProperties
  using Unitful: Pa,hPa, K, °C, km, m, ustrip, uconvert,mbar,atm
  using Unitful
  using Moshi.Match: @match
  abstract type MaterialModel end
  struct Ciddor <: MaterialModel end
  abstract type Mathar <: MaterialModel end
  struct Mathar160_240<: Mathar end
  struct Mathar075_141<: Mathar end
  struct Mathar043_052<: Mathar end
  struct Mathar028_042<: Mathar end
  struct Mathar013_025<: Mathar end


  struct Carlotti <: MaterialModel end
  include("constants.jl")
  include("utilities.jl")
  refractive_index(::G, ::Vararg{T,5}) where {G<:MaterialModel,T<:AbstractFloat} = throw(ArgumentError("MaterialProperties.jl: refractive_index does not support the current model, use $CURRENTMATERIALMODEL instead."))
  refractive_index(model::G,params::Vararg{Int,5}) where {G <: MaterialModel} = refractive_index(model,float.(params)...)
  refractive_index(model::G,params::Vararg{Real,5}) where {G <: MaterialModel} = refractive_index(model,promote.(params)...)
  refractive_index!(::G,::A, ::Vararg{A,5}) where {T<:AbstractFloat,G <: MaterialModel, A<:AbstractArray{T}} = throw(ArgumentError("MaterialProperties.jl: refractive_index! does not support the current model, use $CURRENTMATERIALMODEL instead."))

  # ciddor refractive index model
  include("models/ciddor.jl")
  # mathar refractive index model
  include("models/mathar.jl")
  include("models/carlotti.jl")

  # export models
  export Ciddor,Mathar,Carlotti
  # ciddor model
  export ciddor_refractive_index,ciddor_refractive_index!
  # mathar model
  export mathar013_025_refractive_index,Mathar013_025_refractive_index!
  export mathar028_042_refractive_index,Mathar028_042_refractive_index!
  export mathar043_052_refractive_index,Mathar043_052_refractive_index!
  export mathar075_141_refractive_index,Mathar075_141_refractive_index!
  export mathar160_240_refractive_index,Mathar160_240_refractive_index!
  # carlotti model
  export carlotti_refractive_index,carlotti_refractive_index!
  # refractive index interface general
  export refractive_index,refractive_index!
end
