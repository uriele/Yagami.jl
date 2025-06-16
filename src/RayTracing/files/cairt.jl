using Unitful: ustrip,uconvert,Pa,hPa,km
using NCDatasets: NCDataset

@inline numshort(n::Real;digits=3) = lpad(round(n,digits=digits),digits+4)
@inline numshort(n::Int;digits=3) = lpad(n,digits+4)
@inline textshort(n::String;digits=3) = string(lpad(n,digits+4))

function refractive_info(refractive)
  @info "===================================================================="
  @info "Refractive index information:"
  @info "===================================================================="
  @info "$(textshort("θ_i")) $(textshort("h_j")) $(textshort("(n-1) [1e-4]"))"
  for j in axes(refractive,2)
    for i in axes(refractive,1)
      n= refractive[i,j]
      @info "$(numshort(i)) $(numshort(j)) $(numshort((n-1)*1e4))"
    end
  end
end
function atmosphere_info(atmosphere)
      @info "===================================================================="
      @info "Temperature information:"
      @info "===================================================================="
      @info "$(textshort("θ [°]")) $(textshort("h [km]")) $(textshort("T [K]"))  $(textshort("P [Pa]")) $(textshort("humidity [%]")) $(textshort("CO2 [ppm]")) $(textshort("λ [μm]"))"
      for j in eachindex(atmosphere.temperature.knots_h)
        h= atmosphere.temperature.knots_h[j]
        for i in eachindex(atmosphere.temperature.knots_θ[1:end-1])
          θ= atmosphere.temperature.knots_θ[i]
          T= atmosphere.temperature(θ,h)
          P= atmosphere.pressure(θ,h)
          H= atmosphere.humidity(θ,h)*100
          CO2= atmosphere.co2ppm(θ,h)
          λ= atmosphere.wavelength(θ,h)
          @info "$(numshort(θ)) $(numshort(h)) $(numshort(T)) $(numshort(P)) $(numshort(H)) $(numshort(CO2)) $(numshort(λ))"
        end
      end
end



struct RayTracingProblem{T,MT<:MeanType,AM<:AirModel,EA<:EarthApproximation,ATM<:AtmosphereSetting,V<:AbstractVector{T},M<:AbstractMatrix{T}}
  filename::String
  meantype::MT
  model::AM
  earthmodel::EA
  ###### Data
  atmosphere::ATM
  refractive::M
  pointsx::V
  pointsy::V
  directionsx::V
  directionsy::V
  ##############
  tangent_h:: V
  tangent_θ:: V
  nscans::Int
  nlos::Int

  function RayTracingProblem(filename::String;
    meantype::MT=GeometricMean(), model::AM=Ciddor(),
    earthmodel::EA=Fukushima(),
    logger::L=NullLogger(),
    ) where {MT<:MeanType,AM<:AirModel,
    EA<:EarthApproximation,L<:AbstractLogger}
    nc=NCDataset(filename)
    # number of scans
    nscans = nc.dim["nscans"]
    nlos   = nc.dim["nlos"]

    with_logger(logger) do
      @info "===================================================================="
      @info "Reading data from $filename with model $model and meantype $meantype"
      @info "===================================================================="
      @info "Number of scans: $nscans, Number of lines of sight: $nlos"
      @info " "
    end

    ncscan = nc.group["scans"];
    ncatm  = nc.group["atm"];
    ncorbit= nc.group["orbit"];

    # read the data
    tangent_h = ncscan["tangent_points"][1,:,:]
    tangent_θ = ncscan["tangent_points"][2,:,:]

    # read info
    r_sat = ncorbit["radius"][1]
    θ_sat = ncscan["scans_phi"][:]
    view_angle = ncscan["view_angles"][:,:]

    pointsx = similar(θ_sat,size(view_angle))
    pointsy = similar(pointsx)
    directionsx = similar(pointsx)
    directionsy = similar(pointsy)

    fearth=earthmodelfunction(earthmodel)


    with_logger(logger) do
      @info "Calculating points and directions for the satellite scans"
      @info "$(textshort("los")) $(textshort("scan")) $(textshort("x [km]")) $(textshort("y [km]")) $(textshort("dx")) $(textshort("dy")) $(textshort("view angle [°]")) $(textshort("qt h [km]")) $(textshort("qt θ [°]"))"
    end
    @inbounds for j in axes(pointsx,2)
      pointx =r_sat * cosd(θ_sat[j])
      pointy =r_sat * sind(θ_sat[j])
      h,θ=fearth(pointx,pointy,MAJORAXIS(),MINORAXIS(),COMPLECCENTRICITY²())
      normalx = cosd(θ)
      normaly = sind(θ)
      @inbounds for i in axes(pointsy,1)
        pointsx[i,j] = pointx
        pointsy[i,j] = pointy
        directionsx[i,j],directionsy[i,j] = nadir_angle_normal(normalx,normaly,-view_angle[i,j])
        with_logger(logger) do
          @info "$(numshort(i))  $(numshort(j))  $(numshort(pointx))  $(numshort(pointy))  $(numshort(directionsx[i,j]))  $(numshort(directionsy[i,j]))  $(numshort(view_angle[i,j])) $(numshort(tangent_h[i,j]))  $(numshort(tangent_θ[i,j]))"
        end
      end
    end
    pointsx=pointsx[:]
    pointsy=pointsy[:]
    directionsx=directionsx[:]
    directionsy=directionsy[:]
    tangent_h=tangent_h[:]
    tangent_θ=tangent_θ[:]
    hᵢ = ncatm["z"][:]
    θᵢ = ncatm["phi"][:]

    # not used for now
    hitrancodes =ncatm["hitran"][:]

    pressure = ncatm["pressure"][:,:] |> permutedims |>
      fn-> @. uconvert(Pa, fn*hPa)  |> # convert to hP
      fn-> @. ustrip(fn)
    temperature = ncatm["temperature"][:,:] |> permutedims

    # sort in the reverse order of height
    idxj=similar(hᵢ,Int)
    sortperm!(idxj,hᵢ;rev=true) # sort indices in descending order
    hᵢ= hᵢ[idxj] # sort height by index
    idxi=similar(θᵢ,Int)

    # sort in the range [0,360)
    θᵢ = @. mod(θᵢ, 360) # ensure angles are in [0,360)
    sortperm!(idxi,θᵢ) # sort indices in ascending order
    θᵢ = θᵢ[idxi] # sort angles by index

    # sort the pressure and temperature
    temperature = temperature[idxi,idxj] # sort temperature by index
    pressure = pressure[idxi,idxj] # sort pressure by index

    temperature[ismissing.(temperature)] .= NaN
    pressure[ismissing.(pressure)] .= NaN
    pressure = convert(Matrix{Float64}, pressure) # convert to Float64
    temperature = convert(Matrix{Float64}, temperature) # convert to Float64

    # create the atmosphere object to be used also in CG Integral
    atmosphere=create_atmosphere(;θᵢ=θᵢ,hᵢ=hᵢ,
      temperatureᵢ=temperature,
      pressureᵢ=pressure,
      knots_θ=θᵢ,
      knots_h=hᵢ)


    refractive= grid_refractiveindex(atmosphere;model=Ciddor(),meantype=GeometricMean())
    V=typeof(directionsx)
    M=typeof(refractive)
    T=eltype(directionsx)
    ATM= typeof(atmosphere)

    with_logger(logger) do
      atmosphere_info(atmosphere)
      refractive_info(refractive)
    end


    new{T,MT,AM,EA,ATM,V,M}(
      filename, meantype, model, earthmodel,
      atmosphere,refractive,
      pointsx,pointsy,
      directionsx, directionsy,
      tangent_h,tangent_θ,
      nscans,nlos )
  end
end
Base.show(io::IO,::Type{<:RayTracingProblem{T}}) where T = print(io, "RayTracingProblem{$T}")
Base.show(io::IO, probl::RayTracingProblem{T}) where T = begin
  println("RayTracingProblem{$T}:")
  println("   Filename: $(probl.filename)")
  println("   MeanType: $(typeof(probl.meantype))")
  println("   Model: $(typeof(probl.model))")
  println("   Earth Model: $(typeof(probl.earthmodel))")
  println("   Number of scans: $(probl.nscans)")
  println("   Number of lines of sight: $(probl.nlos)")
end
Base.propertynames(::RayTracingProblem) = PROBLPROPTOTAL

function Base.getproperty(probl::RayTracingProblem, sym::Symbol)
  if sym in PROBLKNOTS
    return getfield(getfield(getfield(probl,:atmosphere),:temperature), sym)
  elseif sym in PROBLATMPROP
    return getfield(getfield(probl,:atmosphere),sym)
  elseif sym in PROBLPROPDIRECT
    return getfield(probl, sym)
  end
  throw(ArgumentError("Property $sym not found in RayTracingProblem."))
end
