


# Interface for Ciddor's refractive index model
refractive_index(::Mathar013_025, temperature::T, pressure::T, wavelength::W=10.0,
humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C} =  mathar013_025_refractive_index(temperature, pressure, wavelength, humidity)
refractive_index!(::Mathar013_025, n::AbstractArray{T}, temperature::AbstractArray{T}, pressure::AbstractArray{T}, wavelength::W=10.0,
humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C} = mathar013_025_refractive_index!(n,temperature, pressure, wavelength, humidity)

refractive_index(::Mathar028_042, temperature::T, pressure::T, wavelength::W=10.0,
humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C} = mathar028_042_refractive_index(temperature, pressure, wavelength, humidity)
refractive_index!(::Mathar028_042, n::AbstractArray{T}, temperature::AbstractArray{T}, pressure::AbstractArray{T},
 wavelength::W=10.0, humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C}=mathar028_042_refractive_index!(n,temperature, pressure, wavelength, humidity)

refractive_index(::Mathar043_052, temperature::T, pressure::T, wavelength::W=10.0,
humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C} =  mathar043_052_refractive_index(temperature, pressure, wavelength, humidity)
refractive_index!(::Mathar043_052, n::AbstractArray{T}, temperature::AbstractArray{T}, pressure::AbstractArray{T}, wavelength::W=10.0,
humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C}=
 mathar043_052_refractive_index!(n,temperature, pressure, wavelength, humidity)

refractive_index(::Mathar075_141, temperature::T, pressure::T, wavelength::W=10.0,
humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C} =  mathar075_141_refractive_index(temperature, pressure, wavelength, humidity)
refractive_index!(::Mathar075_141, n::AbstractArray{T}, temperature::AbstractArray{T}, pressure::AbstractArray{T}, wavelength::W=10.0,
humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C}= mathar075_141_refractive_index!(n,temperature, pressure, wavelength, humidity)

refractive_index(::Mathar160_240, temperature::T, pressure::T, wavelength::W=10.0,
humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C} =  mathar160_240_refractive_index(temperature, pressure, wavelength, humidity)
refractive_index!(::Mathar160_240, n::AbstractArray{T}, temperature::AbstractArray{T}, pressure::AbstractArray{T}, wavelength::W=10.0,
humidity::H=0.0, CO2ppm::C=450.0) where {T<:AbstractFloat,W,H,C}=
 mathar160_240_refractive_index!(n,temperature, pressure, wavelength, humidity)


##################################################################################################################################################
using Unitful


for mathar in MATHARRANGES

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
  REFSIGMA         = Symbol("M",mathar, "REFσ")

  @info "$REFSIGMA"

  @eval function $interfun(temperature::T, pressure::T, wavelength::T, humidity::T) where T<:AbstractFloat
    return $fun(temperature, pressure, wavelength, humidity)
  end
  @eval function $interfun!(n,temperature::AbstractArray{T}, pressure::AbstractArray{T}, wavelength::W, humidity::H) where {T<:AbstractFloat,W,H}
    return $fun!(n,temperature, pressure, wavelength, humidity)
  end

  @eval @inline function $fun(temperature::T, pressure::T, wavelength::T, humidity::T)::T where T<:AbstractFloat
    # Convert temperature to Kelvin if necessary
    σ = CONVERSIONCMTOμM/wavelength # μm^-1 -> cm^-1
    n =one(T)

    # reduced parameters
    TEMPERATURE=(1/temperature-1/MATREF_T)
    HUMIDITY=(humidity*100-MATREF_H) # convert to fraction
    PRESSURE=(pressure-MATREF_P)
    SIGMA = σ - $REFSIGMA

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
    return n
  end

  @eval @inline @generated function $fun!(n::AbstractArray{T},temperature::AbstractArray{T}, pressure::AbstractArray{T},
    wavelength, humidity) where {T<:AbstractFloat}

    expr = Expr[]
    if wavelength <: AbstractArray
      push!(expr, :(wl=wavelength[i]))
    elseif wavelength <: Real
      push!(expr, :(wl=wavelength))
    else
      error("mathar_refractive_index!: wavelength must be a Real or AbstractArray")
    end
    if humidity <: AbstractArray
      push!(expr, :(hu=humidity[i]))
    elseif humidity <: Real
      push!(expr, :(hu=humidity))
    else
      error("mathar_refractive_index!: humidity must be a Real or AbstractArray")
    end

    quote
      for i in eachindex(n)

        $(expr...)

        σ = CONVERSIONCMTOμM/wl # μm^-1 -> cm^-1
        n[i] =one(T)

        # reduced parameters
        TEMPERATURE=(1/temperature[i]-1/MATREF_T)
        HUMIDITY=(hu*100-MATREF_H) # convert to fraction
        PRESSURE=(pressure[i]-MATREF_P)
        SIGMA = σ - $($REFSIGMA)

        TEMPERATURE² = TEMPERATURE*TEMPERATURE
        HUMIDITY² = HUMIDITY*HUMIDITY
        PRESSURE² = PRESSURE*PRESSURE

        TEMPERATURHUMIDITY = TEMPERATURE*HUMIDITY
        TEMPERATUREPRESSURE = TEMPERATURE*PRESSURE
        HUMIDITYPRESSURE = HUMIDITY*PRESSURE
        SIGMAʲ= one(T)

        for j in eachindex($($CP))
          #first step sigmaʲ is 1 (σ^0) then it is multiplied by σ for each loop
          @inbounds n[i] += (
              $($REFC)[j]+
              $($CT)[j]*TEMPERATURE         +$($CTT)[j]*TEMPERATURE²+
              $($CH)[j]*HUMIDITY            +$($CHH)[j]*HUMIDITY²+
              $($CP)[j]*PRESSURE            +$($CPP)[j]*PRESSURE²+
              $($CTH)[j]*TEMPERATURHUMIDITY+
              $($CTP)[j]*TEMPERATUREPRESSURE+
              $($CHP)[j]*HUMIDITYPRESSURE)*SIGMAʲ
          SIGMAʲ *= SIGMA
        end
      end
    end
  end
end
