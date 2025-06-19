using ..YagamiCore: file_to_array

@inline __test_files_existance(inp_folder::String) =begin
  # check all the files necessary exists in the folder
  for (file,header) in (GEOFITFILES,GEOFITHEADERS)
  altitude_file = joinpath(inp_folder, "in_alt.dat")
  if !isfile(altitude_file)
    errormessage = "Altitude file not found: $altitude_file in folder $inp_folder"
    with_logger(logger)
      @error errormessage
      @error "Please ensure that the altitude file is present in the specified folder."
    end
    throw(ArgumentError(errormessage))
  end
end


"""
   $SIGNATURES

Create a `RayTracingProblem` from a Mipas folder.
This struct represents a ray tracing problem with all the necessary data
to perform ray tracing calculations, including the atmosphere, refractive index, and satellite scan information.
# Arguments:
- `filename::String`: The path to the Cairt file containing the ray tracing data.

# Key Arguments:
- `meantype::MT=GeometricMean()`: The mean type to be used for the refractive index calculation (default is `GeometricMean`).
- `model::AM=Ciddor()`: The air model to be used for the refractive index calculation (default is `Ciddor`).
- `earthmodel::EA=Fukushima()`: The earth approximation model to be used (default is `Fukushima`).
- `logger=nothing`: An optional logger to log information during the problem creation (default is `NullLogger`).
# Returns:
- `RayTracingProblem`: An instance of `RayTracingProblem` containing the data read from the Cairt file and the calculated refractive index and atmosphere.
"""
function GeofitRayTracingProblem(folder::String;
    meantype::MT=GeometricMean(), model::AM=Ciddor(),
    earthmodel::EA=Fukushima(),
    logger=nothing,
  )   where {MT<:MeanType,AM<:AirModel,
    EA<:EarthApproximation
  }

  inp_folder=joinpath(folder)

  # check if all necessary files are present in the folder
  __test_files_existance(inp_folder)

  # sanity check for dimensions
  nalt=open(altitude_file, "r") do f
    convert(Int,readlines(f)[22])
  end
  nazimuth=open(azimuth_file, "r") do f
    convert(Int,readlines(f)[16])
  end
  #




  altitude_header = 24
  hᵢ=file_to_array(altitude_file, altitude_header, Float64, " ") # read the altitude file

  azimuth_file = joinpath(inp_folder, "in_lat.dat")
  azimuth_header = 15

  θᵢ=file_to_array(azimuth_file, azimuth_header, Float64, " ") # read the azimuth file

    logger = @match logger begin
      ::AbstractString => infologger(logger)
      ::AbstractLogger => logger
      _ => NullLogger()
    end

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


    hᵢ=earthmodelfunction(earthmodel)


    with_logger(logger) do
      @info "Calculating points and directions for the satellite scans"

      @info "$(textshort("los")) $(textshort("scan")) $(textshort("x [km]")) $(textshort("y [km]")) $(textshort("dx")) $(textshort("dy")) $(textshort("view angle [°]")) $(textshort("qt h [km]")) $(textshort("qt θ [°]"))"
    end
    a² = MAJORAXIS()*MAJORAXIS()
    b² = MINORAXIS()*MINORAXIS()
    @inbounds for j in axes(pointsx,2)
      pointx =r_sat * cosd(θ_sat[j])
      pointy =r_sat * sind(θ_sat[j])
      #h,θ=fearth(pointx,pointy,MAJORAXIS(),MINORAXIS(),COMPLECCENTRICITY²())
      ############################################
      # I do not need theta and since I know that at every point
      # the normal is equal to
      #  n~(x/a^2,y/b^2)/|*|  =>  n~(x/a^2*b^2,y)/|*|
      ############################################
      normalx = pointx/a²
      normaly = pointy/b²
      norm_hypot = hypot(normalx, normaly)
      normalx /= norm_hypot
      normaly /= norm_hypot
      ###########################################
      @inbounds for i in axes(pointsy,1)
        pointsx[i,j] = pointx
        pointsy[i,j] = pointy
        directionsx[i,j],directionsy[i,j] = nadir_angle_normal(-normalx,-normaly,-view_angle[i,j])
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

    L= typeof(logger)

    new{T,MT,AM,EA,ATM,V,M,L}(
      filename, meantype, model, earthmodel,
      atmosphere,refractive,
      pointsx,pointsy,
      directionsx, directionsy,
      tangent_h,tangent_θ,
      nscans,nlos,logger )
end
