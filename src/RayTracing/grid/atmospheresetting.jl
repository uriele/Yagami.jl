"""
   `AtmosphereSetting(knots_θ, knots_h, temperature, pressure, humidity=0.0, co2ppm=0.0, wavelength=10.0)`

Create an atmosphere setting interpolation object for ray tracing.

# Arguments
- `knots_θ::AbstractVector{T}`: A vector of angles in degrees, where `T` is a subtype of `AbstractFloat`.
- `knots_h::AbstractVector{T}`: A vector of heights in kilometers, where `T` is a subtype of `AbstractFloat`.
- `temperature::AbstractMatrix{T}`: A matrix of temperature values, where each row corresponds to a level and each column corresponds to a radius.
- `pressure::AbstractMatrix{T}`: A matrix of pressure values, where each row corresponds to a level and each column corresponds to a radius.
- `humidity::AbstractFloat` (optional): A scalar or matrix of humidity values, default is `0.0`. If a scalar is provided, it is broadcasted to match the size of the temperature matrix.
- `co2ppm::AbstractFloat` (optional): A scalar or matrix of CO2 concentration values in parts per million, default is `0.0`. If a scalar is provided, it is broadcasted to match the size of the temperature matrix.
- `wavelength::AbstractFloat` (optional): A scalar or matrix of wavelength values, default is `10.0`. If a scalar is provided, it is broadcasted to match the size of the temperature matrix.
# Returns
- An `AtmosphereSetting` object containing interpolated atmospheric properties for ray tracing.
"""
struct AtmosphereSetting{N,M,T<:AbstractFloat}
  temperature::AtmInterpolate{N,M,T}
  pressure::AtmInterpolate{N,M,T}
  humidity::AtmInterpolate
  co2ppm::AtmInterpolate
  wavelength::AtmInterpolate
end
@generated function AtmosphereSetting(knots_θ::Vi,knots_h::Vj,
  temperature::M,pressure::M,
  humidity=0.0,co2ppm=0.0,
  wavelength=10.0) where {T<:AbstractFloat,Vi<:AbstractVector{T},Vj<:AbstractVector{T},M<:AbstractMatrix{T}}

  expr=Expr[]

  push!(expr, :(temperature = AtmInterpolate(knots_θ,knots_h,temperature)))
  push!(expr, :(pressure = AtmInterpolate(knots_θ,knots_h,pressure;logh=true)))

  if humidity <: AbstractFloat
    push!(expr, :(humidity = AtmInterpolate([knots_θ[1],knots_θ[end]],[knots_h[1],knots_h[end]],fill(humidity,2,2))))
  elseif humidity <: AbstractMatrix
    push!(expr, :(humidity = AtmInterpolate(knots_θ,knots_h,humidity;logh=true)))
  else
    throw(ArgumentError("humidity must be a AbstractFloat or AbstractMatrix but got $(typeof(humidity))."))
  end

  if co2ppm <: AbstractFloat
    push!(expr, :(co2ppm = AtmInterpolate([knots_θ[1],knots_θ[end]],[knots_h[1],knots_h[end]],fill(co2ppm,2,2))))
  elseif co2ppm <: AbstractMatrix
    push!(expr, :(co2ppm = AtmInterpolate(knots_θ,knots_h,co2ppm;logh=false)))
  else
    throw(ArgumentError("co2ppm must be a AbstractFloat or AbstractMatrix but got $(typeof(co2ppm))."))
  end

  push!(expr, :(wavelength = AtmInterpolate([knots_θ[1],knots_θ[end]],[knots_h[1],knots_h[end]],fill(wavelength,2,2))))

  return quote
    $(expr...)
    return AtmosphereSetting{length(knots_θ),length(knots_h),T}(temperature,pressure,humidity,co2ppm,wavelength)
  end
end

Base.show(io::IO, atm::AtmosphereSetting{N,M,T}) where {N,M,T} = begin
  println(io, "AtmosphereSetting{", N, ",", M, ",", T, "}")

  pressure     = extrema(getfield(atm.pressure,:A))
  temperature  = extrema(getfield(atm.temperature,:A))

  humidity     = extrema(getfield(atm.humidity,:A))   |> unique
  co2ppm       = extrema(getfield(atm.co2ppm,:A))     |> unique
  wavelength   = extrema(getfield(atm.wavelength,:A)) |> unique

  println(io, "  Temperature: ", temperature," K")
  println(io, "  Pressure: ", pressure," Pa")
  println(io, "  Humidity: ", humidity)
  println(io, "  CO2: ", co2ppm, " ppm")
  println(io, "  Wavelength: ", wavelength, " μm")
end

"""
  $SIGNATURES
Create a set of levels for the raytracing method. The levels are evenly spaced in the range `[hmin, hmax]` if `logscale=false`, or logarithmically spaced if `logscale=true`.

## Note
The levels are sorted in descending order, meaning the first element is the maximum level and the last element is the minimum level.
This is the format required by the reay tracing algorithm.

# Arguments
- `T::Type{<:AbstractFloat}`: The type of the levels, default is `Float64`.
- `hmin::T`: The minimum value of the levels, default is `0.0`.
- `hmax::T`: The maximum value of the levels, default is `120.0`.
- `levels::Int`: The number of levels to create, must be a positive integer, default is `100`.
- `logscale::Bool`: If `true`, the levels are logarithmically spaced; if `false`, they are linearly spaced. Default is `true`.
# Returns
- A vector of levels of type `T`, either logarithmically or linearly spaced between `hmin` and `hmax`.
"""
function create_hlevelset(::Type{T}=Float64;hmin::AbstractFloat=0.0,
  hmax::AbstractFloat=120.0,levels::Int=100,logscale=true) where {T<:AbstractFloat}
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
  throw(error("Unexpected error in create_hlevelset."))
end


"""
  $SIGNATURES
Create a set of radii for the raytracing method. The radii are evenly spaced in the range `[θmin, θmax]`.
## Note
The radii are sorted in ascending order, meaning the first element is the minimum radius and the last element is the maximum radius.
# Arguments
- `T::Type{<:AbstractFloat}`: The type of the radii, default is `Float64`.
- `θmin::AbstractFloat`: The minimum angle in degrees, default is `0.0`.
- `θmax::AbstractFloat`: The maximum angle in degrees, default is `359.0`.
- `radii::Int`: The number of radii to create, must be a positive integer, default is `360`.
# Returns
- A vector of radii of type `T`, evenly spaced between `θmin` and `θmax`.
"""
function create_radii(::Type{T}=Float64;θmin::AbstractFloat=0.0,
  θmax::AbstractFloat=359.0,radii::Int=360) where {T<:AbstractFloat}
  @assert(-360<=θmin<=360 ,"θmin must be in the range [-360, 360] but got θmin=$θmin.")
  @assert(-360<=θmax<=360 ,"θmax must be in the range [-360, 360] but got θmax=$θmax.")
  @assert(radii > 0 ,"radii must be a positive integer but got radii=$radii.")
  _θmin = mod(θmin, 360)
  _θmax = mod(θmax, 360)
  @assert (_θmin != _θmax) "θmin must be less than θmax or equal to it, but got θmin=$θmin and θmax=$θmax."
  θ=LinRange(θmin, θmax, radii)
  SVector{radii,T}(θ)
end

"""
 `create_atmosphere(; θᵢ, hᵢ, temperatureᵢ, pressureᵢ, humidity=0.0, co2ppm=0.0, wavelength=10.0, knots_θ=create_hlevelset(T;), knots_h=create_radii(T;))`


Create an atmosphere setting for ray tracing. This function generates an `AtmosphereSetting` object based on the provided parameters, which include angles, heights, temperature, pressure, humidity, CO2 concentration, and wavelength.

# Key Arguments
- `θᵢ::AbstractVector{T}`: A vector of angles in degrees, where `T` is a subtype of `AbstractFloat`.
- `hᵢ::AbstractVector{T}`: A vector of heights in kilometers, where `T` is a subtype of `AbstractFloat`.
- `temperatureᵢ::AbstractMatrix{T}`: A matrix of temperature values, where each row corresponds to a level and each column corresponds to a radius.
- `pressureᵢ::AbstractMatrix{T}`: A matrix of pressure values, where each row corresponds to a level and each column corresponds to a radius.
- `humidity=0.0`: A scalar or matrix of humidity values, default is `0.0`. If a scalar is provided, it is broadcasted to match the size of the temperature matrix.
- `co2ppm=0.0`: A scalar or matrix of CO2 concentration values in parts per million, default is `0.0`. If a scalar is provided, it is broadcasted to match the size of the temperature matrix.
- `wavelength=10.0`: A scalar or matrix of wavelength values, default is `10.0`. If a scalar is provided, it is broadcasted to match the size of the temperature matrix.
- `knots_θ::AbstractVector{T}=create_hlevelset(T;)`: A vector of angles in degrees for the knots, default is created using `create_hlevelset`.
- `knots_h::AbstractVector{T}=create_radii(T;)`: A vector of heights in kilometers for the knots, default is created using `create_radii`.
# Returns
- An `AtmosphereSetting` object containing interpolated atmospheric properties for ray tracing.
"""
@generated function create_atmosphere(;
  θᵢ::AbstractVector{T},
  hᵢ::AbstractVector{T},
  temperatureᵢ::AbstractMatrix{T},
  pressureᵢ::AbstractMatrix{T},
  humidity::T=0.0,
  co2ppm::T=0.0,
  wavelength::T=10.0,
  knots_θ::AbstractVector{T}=create_hlevelset(),
  knots_h::AbstractVector{T}=create_radii(),
) where {T<:AbstractFloat}

  expr = Expr[]
  expr_loop = Expr[]
  # First create interpolation object
  push!(expr,:(atm=AtmosphereSetting(θᵢ,hᵢ,temperatureᵢ,pressureᵢ,humidity,co2ppm, wavelength)))
  push!(expr,:(temperature=similar(temperatureᵢ,length(knots_θ),length(knots_h))))
  push!(expr,:(pressure=similar(pressureᵢ,length(knots_θ),length(knots_h))))


  push!(expr_loop,:(temperature[i,j] = atm.temperature(knots_θ[i],knots_h[j]) ))
  push!(expr_loop,:(pressure[i,j] = atm.pressure(knots_θ[i],knots_h[j]) ))


  if humidity <: AbstractMatrix
    push!(expr,:(humidity = similar(temperatureᵢ,length(knots_θ), length(knots_h))))
    push!(expr_loop,:(humidity[i,j] = atm.humidity(knots_θ[i],knots_h[j]) ))
  elseif humidity <: AbstractFloat
    push!(expr,:(humidity = humidity))
  else
    throw(ArgumentError("humidity must be a AbstractFloat or AbstractMatrix but got $(typeof(humidity))."))
  end

  if co2ppm <: AbstractMatrix
    push!(expr,:(co2ppm = similar(temperatureᵢ,length(knots_θ), length(knots_h))))
    push!(expr_loop,:(co2ppm[i,j] = atm.co2ppm(knots_θ[i],knots_h[j]) ))
  elseif co2ppm <: AbstractFloat
    push!(expr,:(co2ppm = co2ppm))
  else
    throw(ArgumentError("co2ppm must be a AbstractFloat or AbstractMatrix but got $(typeof(co2ppm))."))
  end

  if wavelength <: AbstractMatrix
    push!(expr,:(wavelength = similar(temperatureᵢ,length(knots_θ), length(knots_h))))
    push!(expr_loop,:(wavelength[i,j] = atm.wavelength(knots_θ[i],knots_h[j]) ))
  elseif wavelength <: AbstractFloat
    push!(expr,:(wavelength = wavelength))
  else
    throw(ArgumentError("wavelength must be a AbstractFloat or AbstractMatrix but got $(typeof(wavelength))."))
  end

  return quote
    $(expr...)
    @inbounds for j in eachindex(knots_h)
      @inbounds for i in eachindex(knots_θ)
        $(expr_loop...)
      end
    end
    return AtmosphereSetting(knots_θ,knots_h,temperature,pressure,humidity,co2ppm,wavelength)
  end
end


function grid_refractiveindex(model::AM,mean::MT,atm::AtmosphereSetting{N,M,T}) where {N,M,T<:AbstractFloat,AM,MT}


  n= similar(getfield(atm.temperature,:A),N,M-1)

  knots_θ = getfield(atm.temperature,:knots_θ)
  knots_h = getfield(atm.temperature,:knots_h)

  @inbounds for j in firstindex(knots_h):lastindex(knots_h)-1 # I add extra index for the wrapping
    @inbounds for i in firstindex(knots_θ):lastindex(knots_θ)-1 # I add extra index for the wrapping
      # temperature
      topleft_temp = atm.temperature(knots_θ[i],knots_h[j])
      topright_temp = atm.temperature(knots_θ[i+1],knots_h[j])
      bottomleft_temp = atm.temperature(knots_θ[i],knots_h[j+1])
      bottomright_temp = atm.temperature(knots_θ[i+1],knots_h[j+1])
      # pressure
      topleft_press = atm.pressure(knots_θ[i],knots_h[j])
      topright_press = atm.pressure(knots_θ[i+1],knots_h[j])
      bottomleft_press = atm.pressure(knots_θ[i],knots_h[j+1])
      bottomright_press = atm.pressure(knots_θ[i+1],knots_h[j+1])
      # humidity
      topleft_hum = atm.humidity(knots_θ[i],knots_h[j])
      topright_hum = atm.humidity(knots_θ[i+1],knots_h[j])
      bottomleft_hum = atm.humidity(knots_θ[i],knots_h[j+1])
      bottomright_hum = atm.humidity(knots_θ[i+1],knots_h[j+1])
      # co2ppm
      topleft_co2 = atm.co2ppm(knots_θ[i],knots_h[j])
      topright_co2 = atm.co2ppm(knots_θ[i+1],knots_h[j])
      bottomleft_co2 = atm.co2ppm(knots_θ[i],knots_h[j+1])
      bottomright_co2 = atm.co2ppm(knots_θ[i+1],knots_h[j+1])
      # wavelength
      topleft_wl = atm.wavelength(knots_θ[i],knots_h[j])
      topright_wl = atm.wavelength(knots_θ[i+1],knots_h[j])
      bottomleft_wl = atm.wavelength(knots_θ[i],knots_h[j+1])
      bottomright_wl = atm.wavelength(knots_θ[i+1],knots_h[j+1])
      #
      temp_mean = __mean(mean, topleft_temp, topright_temp, bottomleft_temp, bottomright_temp)
      press_mean = __mean(mean, topleft_press, topright_press, bottomleft_press, bottomright_press)
      hum_mean = __mean(mean, topleft_hum, topright_hum, bottomleft_hum, bottomright_hum)
      co2_mean = __mean(mean, topleft_co2, topright_co2, bottomleft_co2, bottomright_co2)
      wl_mean = __mean(mean, topleft_wl, topright_wl, bottomleft_wl, bottomright_wl)
      # calculate refractive index
      n[i,j] = refractive_index(model,temp_mean, press_mean, wl_mean, hum_mean, co2_mean)
    end
  end

  return n
end


grid_refractiveindex(atm::AtmosphereSetting{N,M,T};model::AM=Ciddor(),meantype::MT=GeometricMean()) where {N,M,T<:AbstractFloat,AM,MT} = grid_refractiveindex(model,meantype, atm)
