
module MaterialProperties
  using DocStringExtensions
  using Unitful: Pa,hPa, K, °C, km, m, ustrip, uconvert,mbar,atm
  using Unitful
  using Moshi.Match: @match
  using StaticArrays: SVector, @SVector
  abstract type AirModel end
  struct Ciddor <: AirModel end
  abstract type Mathar <: AirModel end
  struct Mathar160_240<: Mathar end
  struct Mathar075_141<: Mathar end
  struct Mathar043_052<: Mathar end
  struct Mathar028_042<: Mathar end
  struct Mathar013_025<: Mathar end


  struct Carlotti <: AirModel end
  include("constants.jl")
  include("utilities.jl")


  """
     $SIGNATURES

  Returns the refractive index of air at a given wavelength and temperature, pressure, and humidity.
  The refractive index is calculated using the specified air model and parameters.
  # Arguments
  - `model::G`: The air model to use for the calculation, which must be a subtype of `AirModel`.
  - `temperature::Real`: The temperature in Kelvin.
  - `pressure::Real`: The pressure in Pascals.
  - `wavelength::Real`: The wavelength in micrometers.
  - `humidity::Real`: The relative humidity as a percentage (0-100).
  - `CO2ppm::Real`: The CO2 concentration in parts per million (ppm
  # Returns
  - The refractive index of air at the specified wavelength and conditions.
  """
  refractive_index(model::G,temperature::T, pressure::T, wavelength::T,
  humidity::T, CO2ppm::T) where {G <: AirModel,T<:AbstractFloat} = throw(ArgumentError("MaterialProperties.jl: refractive_index does not support $G, use $CURRENTAirModel instead."))


  """
     $SIGNATURES
  Returns the refractive index of air at a given wavelength and temperature, pressure, and humidity.
  The refractive index is calculated using the specified air model and parameters.
  # Arguments
  - `model::G`: The air model to use for the calculation, which must be a subtype of `AirModel`.
  - `n::A`: An array to store the calculated refractive index values.
  - `temperature::A`: An array of temperatures in Kelvin.
  - `pressure::A`: An array of pressures in Pascals.
  - `wavelength::A`: An array of wavelengths in micrometers.
  - `humidity::A`: An array of relative humidity values as percentages (0-100).
  - `CO2ppm::A`: An array of CO2 concentrations in parts per million

  """
  refractive_index!(model::G,n::A,temperature::A,pressure::A,wavelength::A,humidity::A,CO2ppm::A) where {G <: AirModel,A} = throw(ArgumentError("MaterialProperties.jl: refractive_index! does not support the current model, use $CURRENTAirModel instead."))

  # ciddor refractive index model
  include("models/ciddor.jl")
  # mathar refractive index model
  include("models/mathar.jl")
  include("models/carlotti.jl")

  # export models
  export Ciddor,Mathar,Carlotti,AirModel
  export Mathar160_240,Mathar075_141,Mathar043_052,Mathar028_042,Mathar013_025
  # ciddor model
  export ciddor_refractive_index,ciddor_refractive_index!
  # mathar model
  export mathar013_025_refractive_index,mathar013_025_refractive_index!
  export mathar028_042_refractive_index,mathar028_042_refractive_index!
  export mathar043_052_refractive_index,mathar043_052_refractive_index!
  export mathar075_141_refractive_index,mathar075_141_refractive_index!
  export mathar160_240_refractive_index,mathar160_240_refractive_index!
  # carlotti model
  export carlotti_refractive_index,carlotti_refractive_index!
  # refractive index interface general
  export refractive_index,refractive_index!
end
