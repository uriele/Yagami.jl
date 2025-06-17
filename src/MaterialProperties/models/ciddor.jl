
# Interface for Ciddor's refractive index model
refractive_index!(::Ciddor, n::AbstractArray{T}, temperature::AbstractArray{T}, pressure::AbstractArray{T},
   wavelength::W=10.0, humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C} =
  ciddor_refractive_index!(n, temperature, pressure, wavelength, humidity, CO2ppm)
refractive_index(::Ciddor, temperature::T, pressure::T, wavelength::T=10.0, humidity::T=0.0, CO2ppm::T=450.0) where T<:AbstractFloat =
ciddor_refractive_index(temperature, pressure, wavelength, humidity, CO2ppm)



# Direct call to Ciddor's refractive index model


@generated function ciddor_refractive_index!(n::AbstractArray{T},temperature::AbstractArray{T}, pressure::AbstractArray{T},
  wavelength::W=10.0, humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C}

  expr=Expr[]
  if wavelength<:AbstractArray
    push!(expr,:(wl=wavelength[i]))
  elseif wavelength<:AbstractFloat
    push!(expr,:(wl=wavelength))
  else
    error("ciddor_refractive_index!: wavelength must be a AbstractFloat or AbstractArray")
  end

  if humidity<:AbstractArray
    push!(expr,:(hu=humidity[i]))
  elseif humidity<:AbstractFloat
    push!(expr,:(hu=humidity))
  else
    error("ciddor_refractive_index!: humidity must be a AbstractFloat or AbstractArray")
  end

  if CO2ppm<:AbstractArray
    push!(expr,:(co2=CO2ppm[i]))
  elseif CO2ppm<:AbstractFloat
    push!(expr,:(co2=CO2ppm))
  else
    error("ciddor_refractive_index!: CO2ppm must be a AbstractFloat or AbstractArray")
  end

  return quote
    for i in eachindex(temperature)
      $(expr...) # this will be the wavelength, humidity and CO2ppm for each index
      n[i] = __refractive_index_ciddor(temperature[i], pressure[i], wl, hu, co2)
    end
  end
end

ciddor_refractive_index(temperature::T, pressure::T, wavelength::W=10.0, humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C}  = __refractive_index_ciddor(temperature, pressure, wavelength, humidity, CO2ppm)



@inline function __compressibility_ciddor(temperature_celsius::AbstractFloat, temperature_kelvin::AbstractFloat, pressure::AbstractFloat, xw::AbstractFloat)
  fraction = pressure/temperature_kelvin
  fraction² = fraction*fraction
  xw² = xw*xw
  1-fraction*(
    CID_COMPRESSIBILITYA0+
    CID_COMPRESSIBILITYA1*temperature_celsius+
    CID_COMPRESSIBILITYA2*temperature_celsius^2+
    (CID_COMPRESSIBILITYB0+CID_COMPRESSIBILITYB1*temperature_celsius)*xw +
    (CID_COMPRESSIBILITYC0+CID_COMPRESSIBILITYC1*temperature_celsius) *xw²
  )+
  fraction²*(
    CID_COMPRESSIBILITYD+
    CID_COMPRESSIBILITYE*xw²
  )
end

@inline function __density_ciddor(temperature_kelvin::T,pressure::T,molar_mass::T,compressibility::T)::T where T<:AbstractFloat
  pressure* molar_mass/(RGAS*temperature_kelvin*compressibility)
end


@inline function __saturation_vapor_pressure(temperature_kelvin::T)::T where T
  if temperature_kelvin >= CONVERSIONCELSIUSTOKELVIN
    exp(CID_A*temperature_kelvin*temperature_kelvin+CID_B*temperature_kelvin+CID_C+CID_D/temperature_kelvin)
  else
    10^(CID_E/temperature_kelvin+CID_F)
  end
end


##########################################################################################################################
#M_dry_air= MOLAR_MASS_STANDARD_AIR+MOLAR_MASS_CO2_CONTRIBUTION*(CO2ppm-MOLAR_MASS_STANDARD_CO2PPM) # molar mass of dry air, g/mol
const CIDREF_ZAIR  ::Float64 = __compressibility_ciddor(CIDREF_TEMPERATURECELSIUSAIR,CIDREF_TEMPERATUREKELVINAIR,CIDREF_PRESSUREAIR,0.0) # 0% humidity
#M_water = T(18.01528e-3) # molar mass of water vapor, g/mol
const CIDREF_ZWA ::Float64 = __compressibility_ciddor(CIDREF_TEMPERATUREWATERCELSIUS,CIDREF_TEMPERATUREWATERKELVIN,CIDREF_PRESSUREWATER,1.0) # 100% humidity

const CIDREF_RHOWV = __density_ciddor(CIDREF_TEMPERATUREWATERKELVIN, CIDREF_PRESSUREWATER, CIDREF_MOLARMASSW, CIDREF_ZWA) # density of water vapor at standard conditions, g/m^3

@inline function __refractive_index_ciddor(temperature_kelvin::T,pressure::T,wavelength::T,humidity::T,co2ppm::T)::T where T<:AbstractFloat

  temperature_celsius = kelvin_to_celsius(temperature_kelvin)::T # convert temperature to Celsius
  σ² = T(1)/(wavelength*wavelength) # vacuum wavenumber µm^-2
  svp= __saturation_vapor_pressure(temperature_kelvin) # saturation vapor pressure in Pa

  # enhancment factor of water vapor in air
  f = CID_α +CID_β*pressure + CID_γ*temperature_celsius^2

  # molar factor of water vapor in air
  xw = f*humidity*svp/pressure # water vapor molar fraction

  # refractive index of standard air at 15°C and 101325 Pa 0% humidity and 450 ppm CO2
  nas = 1+(CID_K1/(CID_K0-σ²) + CID_K3/(CID_K2-σ²))*1e-8  #standard dry air (n_as-1)

  # refractive index of standard dry air at 15°C and 101325 Pa 0% humidity and co2ppm ppm CO2
  naxs= 1+(nas-1)*(1+CIDREF_CO2*(co2ppm-CIDREF_CO2PPM))

  # refractive index of standard water vapor at 20°C and 1333 Pa 100% humidity
  nws = 1+CID_WVCORRECTING*(CID_W0+σ²*(CID_W1 + σ²*(CID_W2 + CID_W3*σ²))) #

  # molar mass of dry air, kg/mol
  Ma = T(1e-3)*(CIDREF_MOLARMASSAIR+CIDREF_MOLARMASSCO2*(co2ppm-CIDREF_MOLARMASSC02PPM)) # only one that is not a constant but depends on co2ppm
  # molar mass of water vapor, kg/mol
  Mw = CIDREF_MOLARMASSW # it is a constant but I put it here for clarity
  Za = CIDREF_ZAIR # compressibility of dry air at standard conditions (15°C and 101325 Pa) 0% humidity
  Zw = CIDREF_ZWA # compressibility of water vapor at standard conditions (20°C and 1333 Pa) 100% humidity

  # Eq. 4 with (T,P,xw) = (288.15K,101325Pa,0)
  ρaxs = __density_ciddor(CIDREF_TEMPERATUREKELVINAIR, CIDREF_PRESSUREAIR, Ma, Za) # density of dry air at standard conditions, g/m^3

   # Eq. 4 with (T,P,xw) = (293.15K,1333Pa,1)
  ρws = CIDREF_RHOWV

  Zeff = __compressibility_ciddor(temperature_celsius,temperature_kelvin,pressure,xw) # effective compressibility

  ρa = __density_ciddor(temperature_kelvin, pressure, Ma, Zeff)*(one(T)-xw) # density of dry air, g/m^3
  ρw = __density_ciddor(temperature_kelvin, pressure, Mw, Zeff)*xw # density of water vapor, g/m^3

  return (1+(ρa/ρaxs)*(naxs-1)+(ρw/ρws)*(nws-1))::T

end
