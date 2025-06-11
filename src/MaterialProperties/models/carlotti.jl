

# Interface for Ciddor's refractive index model
refractive_index!(::Carlotti, n::AbstractArray{T}, temperature::AbstractArray{T}, pressure::AbstractArray{T}, wavelength::AbstractArray{T},
humidity::AbstractArray{T}, CO2ppm::AbstractArray{T}) where {T<:AbstractFloat} =
  carlotti_refractive_index!(n, temperature, pressure, wavelength, humidity, CO2ppm)
refractive_index(::Carlotti, temperature::T, pressure::T, wavelength::T=10.0, humidity::T=0.0, CO2ppm::T=450.) where T<:AbstractFloat =
  carlotti_refractive_index(temperature, pressure, wavelength, humidity, CO2ppm)



# Direct call to Carltti's refractive index model
function carlotti_refractive_index!(n::A,temperature::A, pressure::A, wavelength::Any=nothing, humidity::Any=0.0, CO2ppm::Any=450.) where {T<:AbstractFloat,A<:AbstractArray{T}}
  @simd for i in eachindex(n)
    @inbounds n[i] = __refractive_index_carlotti(temperature[i], pressure[i])
  end
  return n
end

carlotti_refractive_index(temperature::T, pressure::T, ::Any=nothing, wavelength::Any=0.0, CO2ppm::Any=450.) where T<:AbstractFloat  = __refractive_index_carlotti(temperature, pressure)





@inline function __refractive_index_carlotti(temperature_kelvin::T, pressure_pascal::T) where T<:AbstractFloat
  return T(1.0)+CARLOTTI_BASE*(pressure_pascal/temperature_kelvin)*CARLOTTIREF_PRESSTEMP
end
