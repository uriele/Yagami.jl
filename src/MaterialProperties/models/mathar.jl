


# Interface for Ciddor's refractive index model
refractive_index!(::Mathar, n::A, temperature::A, pressure::A, wavelength::A, humidity::A, CO2ppm::Union{A,Nothing}=450.) where {T<:AbstractFloat,A<:AbstractArray{T}} =
  mathar_refractive_index!(n, temperature, pressure, wavelength, humidity, CO2ppm)
refractive_index(::Mathar, temperature::T, pressure::T, wavelength::T, humidity::T, CO2ppm::Any=nothing) where T<:AbstractFloat =
  mathar_refractive_index(temperature, pressure, wavelength, humidity, CO2ppm)

refractive_index(::Mathar013_025, temperature::T, pressure::T, wavelength::T, humidity::T, CO2ppm::Any=450.) where T<:AbstractFloat =  Mathar013_025_refractive_index(temperature, pressure, wavelength, humidity, CO2ppm)
refractive_index!(::Mathar013_025, n::A, temperature::A, pressure::A, wavelength::A, humidity::A, CO2ppm::Union{A,Nothing}=nothing) where {T<:AbstractFloat,A<:AbstractArray{T}}=  Mathar013_025_refractive_index!(temperature, pressure, wavelength, humidity, CO2ppm)

refractive_index(::Mathar028_042, temperature::T, pressure::T, wavelength::T, humidity::T, CO2ppm::Any=450.) where T<:AbstractFloat =  Mathar028_042_refractive_index(temperature, pressure, wavelength, humidity, CO2ppm)
refractive_index!(::Mathar028_042, n::A, temperature::A, pressure::A, wavelength::A, humidity::A, CO2ppm::Union{A,Nothing}=nothing) where {T<:AbstractFloat,A<:AbstractArray{T}}=  Mathar028_042_refractive_index!(temperature, pressure, wavelength, humidity, CO2ppm)

refractive_index(::Mathar043_052, temperature::T, pressure::T, wavelength::T, humidity::T, CO2ppm::Any=450.) where T<:AbstractFloat =  Mathar043_052_refractive_index(temperature, pressure, wavelength, humidity, CO2ppm)
refractive_index!(::Mathar043_052, n::A, temperature::A, pressure::A, wavelength::A, humidity::A, CO2ppm::Union{A,Nothing}=nothing) where {T<:AbstractFloat,A<:AbstractArray{T}}=  Mathar043_052_refractive_index!(temperature, pressure, wavelength, humidity, CO2ppm)

refractive_index(::Mathar075_141, temperature::T, pressure::T, wavelength::T, humidity::T, CO2ppm::Any=450.) where T<:AbstractFloat =  Mathar075_141_refractive_index(temperature, pressure, wavelength, humidity, CO2ppm)
refractive_index!(::Mathar075_141, n::A, temperature::A, pressure::A, wavelength::A, humidity::A, CO2ppm::Union{A,Nothing}=nothing) where {T<:AbstractFloat,A<:AbstractArray{T}}=  Mathar075_141_refractive_index!(temperature, pressure, wavelength, humidity, CO2ppm)

refractive_index(::Mathar160_240, temperature::T, pressure::T, wavelength::T, humidity::T, CO2ppm::Any=450.) where T<:AbstractFloat =  Mathar160_240_refractive_index(temperature, pressure, wavelength, humidity, CO2ppm)
refractive_index!(::Mathar160_240, n::A, temperature::A, pressure::A, wavelength::A, humidity::A, CO2ppm::Any=450.) where {T<:AbstractFloat,A<:AbstractArray{T}}=  Mathar160_240_refractive_index!(temperature, pressure, wavelength, humidity, CO2ppm)

##################################################################################################################################################
using Unitful


for mathar in ("013_025","028_042","043_052","075_141","160_240")

  modelname         = Symbol("Mathar", mathar)
  interfun         = Symbol("mathar", mathar, "_refractive_index")
  interfun!        = Symbol(interfun, "!")
  lowercase_mathar = lowercase(mathar)
  fun              = Symbol("__refractive_index_mathar", lowercase_mathar)
  fun!             = Symbol(fun, "!")
  REFC             = Symbol("M",mathar, "REFC")
  CT               = Symbol("M",mathar, "CT")
  CTT              = Symbol("M",mathar, "CTT")
  CH               = Symbol("M",mathar, "CH")
  CHH              = Symbol("M",mathar, "CHH")
  CP               = Symbol("M",mathar, "CP")
  CPP              = Symbol("M",mathar, "CPP")
  CTH              = Symbol("M",mathar, "CTH")
  CTP              = Symbol("M",mathar, "CTP")
  CHP              = Symbol("M",mathar, "CHP")
  REFσ             = Symbol("M",mathar, "REFσ")

  @eval function $interfun(temperature::T, pressure::T, wavelength::T, humidity::T, CO2ppm::Any=450.) where T<:AbstractFloat
    return $fun(temperature, pressure, wavelength, humidity)
  end

  @eval @inline function $fun(temperature::T, pressure::T, wavelength::T, humidity::T)::T where T<:AbstractFloat
    # Convert temperature to Kelvin if necessary
    σ = CONVERSIONCMTOμM/wavelength # μm^-1 -> cm^-1
    n =one(T)

    # reduced parameters
    TEMPERATURE=(1/temperature-1/MATREF_T)
    HUMIDITY=(humidity*100-MATREF_H) # convert to fraction
    PRESSURE=(pressure-MATREF_P)
    SIGMA = σ - $REFσ

    TEMPERATURE² = TEMPERATURE*TEMPERATURE
    HUMIDITY² = HUMIDITY*HUMIDITY
    PRESSURE² = PRESSURE*PRESSURE

    TEMPERATURHUMIDITY = TEMPERATURE*HUMIDITY
    TEMPERATUREPRESSURE = TEMPERATURE*PRESSURE
    HUMIDITYPRESSURE = HUMIDITY*PRESSURE
    SIGMAʲ= one(T)

    @inbounds for j in eachindex($CP)
      #first step sigmaʲ is 1 (σ^0) then it is multiplied by σ for each loop
      n += ($REFC[j]+$CT[j]*TEMPERATURE+$CTT[j]*TEMPERATURE²+
            $CH[j]*HUMIDITY+$CHH[j]*HUMIDITY²+
            $CP[j]*PRESSURE+$CPP[j]*PRESSURE²+
            $CTH[j]*TEMPERATURHUMIDITY+
            $CTP[j]*TEMPERATUREPRESSURE+
            $CHP[j]*HUMIDITYPRESSURE)*SIGMAʲ
      SIGMAʲ *= SIGMA
    end
  end

  @eval function $fun!(n::A,temperature::A, pressure::A, wavelength::A, humidity::A)::A where {T<:AbstractFloat,A<:AbstractArray{T}}
    for i in eachindex(n)
      σ = CONVERSIONCMTOμM/wavelength # μm^-1 -> cm^-1
      n[i] =one(T)

      # reduced parameters
      TEMPERATURE=(1/temperature-1/MATREF_T)
      HUMIDITY=(humidity*100-MATREF_H) # convert to fraction
      PRESSURE=(pressure-MATREF_P)
      SIGMA = σ - $REFσ

      TEMPERATURE² = TEMPERATURE*TEMPERATURE
      HUMIDITY² = HUMIDITY*HUMIDITY
      PRESSURE² = PRESSURE*PRESSURE

      TEMPERATURHUMIDITY = TEMPERATURE*HUMIDITY
      TEMPERATUREPRESSURE = TEMPERATURE*PRESSURE
      HUMIDITYPRESSURE = HUMIDITY*PRESSURE
      SIGMAʲ= one(T)

      for j in eachindex($CP)
        #first step sigmaʲ is 1 (σ^0) then it is multiplied by σ for each loop
        @inbounds n[i] += ($REFC[j]+$CT[j]*TEMPERATURE+$CTT[j]*TEMPERATURE²+
            $CH[j]*HUMIDITY+$CHH[j]*HUMIDITY²+
            $CP[j]*PRESSURE+$CPP[j]*PRESSURE²+
            $CTH[j]*TEMPERATURHUMIDITY+
            $CTP[j]*TEMPERATUREPRESSURE+
            $CHP[j]*HUMIDITYPRESSURE)*SIGMAʲ
        SIGMAʲ *= SIGMA
      end

    end
  end
end
