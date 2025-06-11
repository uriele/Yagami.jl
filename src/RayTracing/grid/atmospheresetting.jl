import NaturalNeighbours.interpolate as NNinterpolate
const SMALLSHIFT= 1e-6


"""
  $SIGNATURES
Create a set of levels for the raytracing method. The levels are evenly spaced in the range `[hmin, hmax]` if `logscale=false`, or logarithmically spaced if `logscale=true`.

## Note
The levels are sorted in descending order, meaning the first element is the maximum level and the last element is the minimum level.
This is the format required by the reay tracing algorithm.

# Arguments
- `T::Type{<:Real}`: The type of the levels, default is `Float64`.
- `hmin::T`: The minimum value of the levels, default is `0.0`.
- `hmax::T`: The maximum value of the levels, default is `120.0`.
- `levels::Int`: The number of levels to create, must be a positive integer, default is `100`.
- `logscale::Bool`: If `true`, the levels are logarithmically spaced; if `false`, they are linearly spaced. Default is `true`.
# Returns
- A vector of levels of type `T`, either logarithmically or linearly spaced between `hmin` and `hmax`.
"""
function create_hlevelset(::Type{T}=Float64;hmin::Real=T(0.0),
  hmax::Real=T(120.0),levels::Int=100,logscale=true) where {T<:Real}
  @assert(levels > 0 ,"levels must be a positive integer but got levels=$levels.")
  @assert(hmin >= 0   ,"hmin must be non-negative but got hmin=$hmin.")
  @assert(hmax>hmin  ,"hmax must be greater than hmin but got hmin=$hmin and hmax=$hmax.")
  if logscale
    # shift to avoid log(0)
    shift =zero(T)
    if hmin < SMALLSHIFT
      hmin = SMALLSHIFT
      hmax += SMALLSHIFT
      shift = SMALLSHIFT
    end
    h = Base.LogRange(hmax, hmin, levels).-shift
    return SVector{levels,T}(h)
  else
    h = LinRange(hmax,hmin, levels)
    return SVector{levels,T}(h)
  end

end


"""
  $SIGNATURES
Create a set of radii for the raytracing method. The radii are evenly spaced in the range `[θmin, θmax]`.
## Note
The radii are sorted in ascending order, meaning the first element is the minimum radius and the last element is the maximum radius.
# Arguments
- `T::Type{<:Real}`: The type of the radii, default is `Float64`.
- `θmin::Real`: The minimum angle in degrees, default is `0.0`.
- `θmax::Real`: The maximum angle in degrees, default is `359.0`.
- `radii::Int`: The number of radii to create, must be a positive integer, default is `360`.
# Returns
- A vector of radii of type `T`, evenly spaced between `θmin` and `θmax`.
"""
function create_radii(::Type{T}=Float64;θmin::Real=T(0.0),
  θmax::Real=T(359.0),radii::Int=360) where {T<:Real}
  @assert(-360<=θmin<=360 ,"θmin must be in the range [-360, 360] but got θmin=$θmin.")
  @assert(-360<=θmax<=360 ,"θmax must be in the range [-360, 360] but got θmax=$θmax.")
  @assert(radii > 0 ,"radii must be a positive integer but got radii=$radii.")
  _θmin = mod(θmin, 360)
  _θmax = mod(θmax, 360)
  @assert (_θmin != _θmax) "θmin must be less than θmax or equal to it, but got θmin=$θmin and θmax=$θmax."
  h=LinRange(θmin, θmax, radii)
  SVector{radii,T}(h)
end

"""
  $SIGNATURES
Create an atmosphere model for ray tracing. The atmosphere is defined by the levels, radii, temperature, pressure, humidity, CO2 concentration, and wavelength.
# Arguments
- `levels::AL`: A vector of levels, where `AL` is an abstract vector type.
- `radii::AR`: A vector of radii, where `AR` is an abstract vector type.
- `temperature::MA`: A matrix of temperature values, where `MA` is an abstract matrix type.
- `pressure::MA`: A matrix of pressure values, where `MA` is an abstract matrix type.
- `humidity::MA`: A matrix of humidity values, where `MA` is an abstract matrix type.
- `co2ppm::MA`: A matrix of CO2 concentration values in parts per million, where `MA` is an abstract matrix type.
- `wavelength::MA`: A matrix of wavelength values, where `MA` is an abstract matrix type.
- `model::MO`: An air model type, where `MO` is a subtype of `AirModel`. Default is `Ciddor()`.
# Returns
- An atmosphere model represented as a matrix of refractive indices, where each element corresponds to the refractive index at a specific level and radius.
"""
function create_atmosphere(levels::AL,
  radii::AR, temperature::MA,
  pressure::MA,humidity::MA,co2ppm::MA, wavelength::MA,
  model::MO=Ciddor()) where {T<:Real,MA<:AbstractMatrix{T},AL<:AbstractVector{T},AR<:AbstractVector{T},MO<:AirModel}
  ###############################################################################################
  @assert(length(levels) > 0, "levels must be a non-empty vector but got levels=$levels.")
  @assert(length(radii) > 0, "radii must be a non-empty vector but got radii=$radii.")
  Nt,Mt = size(temperature)
  Np,Mp = size(pressure)
  @assert(Nt == Np && Mt == Mp, "temperature and pressure must have the same dimensions but got temperature=$temperature and pressure=$pressure.")
  N,M=Nt,Mt
  @assert(length(levels) == M+1, "levels must have length M+1 where M is the number of radii but got levels=$(length(levels)) and M=$M.")
  @assert(length(radii) == N, "radii must have length N where N is the number of levels but got radii=$(length(radii)) and N=$N.")
  ###############################################################################################
  atmosphere=similar(temperature)

  refractive_index!(model,atmosphere, temperature, pressure, wavelength, humidity,co2ppm)
  # Add CO2 effect if needed
  return atmosphere

end



"""
  $SIGNATURES
Create an atmosphere model for ray tracing. The atmosphere is defined by the levels, radii, temperature, pressure, humidity, CO2 concentration, and wavelength.
# KeyArguments
- `levels::AbstractVector`: A vector of levels, where each element represents a height level in the atmosphere.
- `radii::AbstractVector`: A vector of radii, where each element represents a radius in the atmosphere.
- `temperature::AbstractMatrix`: A matrix of temperature values, where each row corresponds to a level and each column corresponds to a radius.
- `pressure::AbstractMatrix`: A matrix of pressure values, where each row corresponds to a level and each column corresponds to a radius.
- `hunidity::T` (optional): A scalar or matrix of humidity values, default is `0.0`. If a scalar is provided, it is broadcasted to match the size of the temperature matrix.
- `co2ppm::T` (optional): A scalar or matrix of CO2 concentration values in parts per million, default is `0.0`. If a scalar is provided, it is broadcasted to match the size of the temperature matrix.
- `wavelength::T` (optional): A scalar or matrix of wavelength values, default is `10.0`. If a scalar is provided, it is broadcasted to match the size of the temperature matrix.
- `model::MO` (optional): An air model type, where `MO` is a subtype of `AirModel`. Default is `Ciddor()`.
# Returns
- An atmosphere model represented as a matrix of refractive indices, where each element corresponds to the refractive index at a specific level and radius.
"""
function create_atmosphere(::Type{T}=Float64;
  levels::AbstractVector,
  radii::AbstractVector,
  temperature::AbstractMatrix,
  pressure::AbstractMatrix,
  humidity=T(0.0),
  co2ppm=T(0.0),
  wavelength=T(10.0),
  model::MO=Ciddor()) where {T<:Real,MO<:AirModel}

  levels = SVector{length(levels),T}(levels)
  radii = SVector{length(radii),T}(radii)
  temperature = Matrix{T}(temperature)
  pressure = Matrix{T}(pressure)
  if isa(humidity,Real)
    humidity = fill(humidity, size(temperature))
  else
    humidity = Matrix{T}(humidity)
  end
  if isa(co2ppm,Real)
    co2ppm = fill(co2ppm, size(temperature))
  else
    co2ppm = Matrix{T}(co2ppm)
  end
  if isa(wavelength,Real)
    wavelength = fill(wavelength, size(temperature))
  else
    wavelength = Matrix{T}(wavelength)
  end

  return create_atmosphere(levels, radii, temperature, pressure, humidity, co2ppm, wavelength, model)
end




function interpolate_atmosphere(temperature_values::V,pressure_values::V,
  level_values::V,radius_values::V,
  level_knots::LK, radius_knots::RK,model::G=Ciddor()) where {T<:Real,V<:AbstractVector{T}, LK<:AbstractVector{T}, RK<:AbstractVector{T}, G<:AirModel}




end


function interpolate_atmosphere(temperature_values::M,pressure_values::M,
  level_values::LV,radius_values::RV,
  level_knots::LK, radius_knots::RK,model::G=Ciddor()) where {T<:Real,M<:AbstractMatrix{T},
  LV<:AbstractVector{T}, RV<:AbstractVector{T},
  LK<:AbstractVector{T}, RK<:AbstractVector{T}, G<:AirModel}

  atmosphere = Matrix{T}(undef, length(LK), length(RK)-1)







end



function interpolate_atmosphere(temperature_values::M,pressure_values::M,
  level_values::M,radius_values::M,
  levels::Int,radii::Int;
  logscale=true,
  hmin=zero(T),hmax=T(120),
  θmin=zero(T),θmax=T(359),
  model::MO=Ciddor()) where {T<:Real,MO<:AirModel, M<:AbstractMatrix{T}}
end
