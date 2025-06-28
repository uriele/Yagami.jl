using ..YagamiCore: file_to_array,vmr_to_namedtuple
using Logging,LoggingExtras
using Moshi.Match:@match
using DocStringExtensions
@inline function __test_files_existance(inp_folder::S,logger::L=NullLogger()) where {S<:AbstractString,L<:AbstractLogger}

  # check all the files necessary exists in the folder
  for (file) in (GEOFITFILES)
    input_file = joinpath(inp_folder, file)
    if !isfile(input_file)
      errormessage = "file not found: $file in folder $inp_folder"
      with_logger(logger) do
        @error errormessage
        @error "Please ensure that the altitude file is present in the specified folder."
      end
      throw(ArgumentError(errormessage))
    end
  end
end


@inline function __read_orbit_file(filename::String, logger::L=NullLogger()) where {L<:AbstractLogger}
  # read the orbit file and return the radius and azimuth angles
  open(filename, "r") do io
    lines = readlines(io)
    with_logger(logger) do
      @info "Reading orbit file: $filename"

      str= textshort("X [km]",padding=15) *
           textshort("Y [km]",padding=15) *
           textshort("Limb ang [°]",padding=15) *
           textshort("dirX",padding=15) *
           textshort("dirY",padding=15) *
           textshort("tq_h [km]",padding=15) *
           textshort("tq_θ [°]",padding=15)
      @info str
    end

    # find all sequence lines
    pos=findall(x->!isnothing(match(r"^\sSEQUENCE",x)),lines);

    nscans = length(pos)
    # assume same number of lines of sight for each scan
    nlos   = parse(Int,lines[pos[1]+6])

    orbit = zeros(Float64, 6, nlos* nscans)

    @inbounds for i in eachindex(pos)
       @inbounds for j in 1:nlos
        sj=pos[i]+7+j
        input_read=parse.(Float64,split(lines[sj]))

        pointx = input_read[2] # x coordinate of the point
        pointy = input_read[1] # y coordinate of the point
        angle = input_read[3] # limb angle in degrees
        dirx,diry=limb_angle_normal(-pointx, -pointy, -angle) # calculate the limb direction normal

        with_logger(logger) do
          str = numshortf(pointx,padding=10,digits=8) *
                numshortf(pointy,padding=10,digits=8) *
                numshortf(angle,padding=10,digits=8) *
                numshortf(dirx,padding=10,digits=8) *
                numshortf(diry,padding=10,digits=8) *
                numshortf(input_read[4],padding=15) *
                numshortf(input_read[5],padding=15)
          @info str
        end

        orbit[:,(i-1)*nlos+j] = [pointx,pointy, dirx, diry, input_read[4], input_read[5]]
      end
    end
    orbit,nscans,nlos
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
    humidity::Real=0.0, # used to override the humidity in the file
    co2ppm::Real=400.0,
  )   where {MT<:MeanType,AM<:AirModel,
    EA<:EarthApproximation
  }

  inp_folder=joinpath(folder,"INP_FILES")


  logger = @match logger begin
      ::AbstractString => infologger(logger)
      ::AbstractLogger => logger
      _ => NullLogger()
  end


  # check if all necessary files are present in the folder
  __test_files_existance(inp_folder,logger)



  filename= joinpath(inp_folder,GEOFITFILES[1])
  hᵢ = file_to_array(filename, GEOFITHEADERS[1]) # read the input file
  filename= joinpath(inp_folder,GEOFITFILES[2])
  θᵢ = file_to_array(filename, GEOFITHEADERS[2]) # read the input file
  filename= joinpath(inp_folder,GEOFITFILES[3])
  pressureᵢ = file_to_array(joinpath(inp_folder, GEOFITFILES[3]), GEOFITHEADERS[3]) # read the input file
  filename= joinpath(inp_folder,GEOFITFILES[4])
  temperatureᵢ = file_to_array(joinpath(inp_folder, GEOFITFILES[4]), GEOFITHEADERS[4]) # read the input file


  naltitude = length(hᵢ)
  nazimuth = length(θᵢ)

  temperatureᵢ = reshape(temperatureᵢ, naltitude, nazimuth) |> permutedims
  pressureᵢ = reshape(pressureᵢ, naltitude, nazimuth) |> permutedims |>
      fn-> @. uconvert(Pa, fn*hPa)  |> # convert to Pa
      fn-> @. ustrip(fn)

  #@eval $(GEOFITVARNAMES[end])=vrm_to_namedtuple(joinpath(inp_folder, GEOFITFILES[end]), GEOFITHEADERS[end])
  knots_h = file_to_array(joinpath(inp_folder, "in_levels.dat"), 0)
  knots_θ = file_to_array(joinpath(inp_folder, "in_radii.dat"), 0)
  sort!(knots_h;rev=true) # sort the knots in descending order
  knots_θ.-= 90 # convert the azimuth angles to the range [-90, 90]
  knots_θ= mod.(knots_θ, 360) # ensure the azimuth angles are in the range [0, 360)
  sort!(knots_θ) # sort the azimuth angles in descending order

  idx_theta = similar(θᵢ,Int)
  θᵢ .-=90
  θᵢ   = mod.(θᵢ, 360) # ensure the azimuth angles are in the range [0, 360)
  sortperm!(idx_theta, θᵢ) # sort the azimuth angles in descending order
  θᵢ = θᵢ[idx_theta] # sort the azimuth angles in descending order
  idx_h = similar(hᵢ,Int)
  sortperm!(idx_h, hᵢ;rev=true) # sort the altitude points in descending order
  hᵢ = hᵢ[idx_h] # sort the altitude points in descending order
  temperatureᵢ = temperatureᵢ[idx_theta,idx_h] # sort the temperature profile according to the azimuth angles
  pressureᵢ = pressureᵢ[idx_theta,idx_h] # sort the pressure profile according to the azimuth angles

  with_logger(logger) do
    @info "Read $(length(hᵢ)) altitude points and $(length(θᵢ)) azimuth points."
    @info "Temperature profile: $(size(temperatureᵢ))"
    @info "Pressure profile: $(size(pressureᵢ))"
  end

  with_logger(logger) do
    @info "Read $(length(hᵢ)) altitude points and $(length(θᵢ)) azimuth points."
    @info "Temperature profile: $(size(temperatureᵢ))"
    @info "Pressure profile: $(size(pressureᵢ))"
    for (h, θ, temperature, pressure) in zip(hᵢ, θᵢ, temperatureᵢ, pressureᵢ)
      str= numshort(h,padding=10) *
           numshort(θ,padding=10) *
           numshortf(temperature,padding=10) *
           numshortf(pressure,padding=10)
      @info str
    end

    @info "Knots for altitude: $(knots_h)"
    @info "Knots for azimuth: $(knots_θ)"
    for (h, θ) in zip(knots_h, knots_θ)
      str = numshort(h,padding=10) * numshort(θ,padding=10)
      @info str
    end

  end



  atmosphere=create_atmosphere(;θᵢ=θᵢ,hᵢ=hᵢ,
    temperatureᵢ=temperatureᵢ,
    pressureᵢ=pressureᵢ,
    knots_θ=knots_θ,
    knots_h=knots_h)
  refractive= grid_refractiveindex(atmosphere;model=Ciddor(),meantype=GeometricMean())

  #####################################################################
  orbit,nscans,nlos=__read_orbit_file(joinpath(inp_folder, "orbit.dat"), logger)
  # orbit file points
  pointsx=orbit[1,:] # x coordinates of the points
  pointsy=orbit[2,:] # y coordinates of the points
  directionsx=orbit[3,:] # x components of the directions
  directionsy=orbit[4,:] # y components of the directions
  tangent_h=orbit[5,:] # tangent heights
  tangent_θ=orbit[6,:] # tangent azimuths
  #####################################################################
  V=typeof(directionsx)
  M=typeof(refractive)
  T=eltype(directionsx)
  ATM= typeof(atmosphere)

  with_logger(logger) do
    atmosphere_info(atmosphere)
    refractive_info(refractive)
  end


  L= typeof(logger)


  RayTracingProblem{T,MT,AM,EA,ATM,V,M,L}(
      filename, meantype, model, earthmodel,
      atmosphere,refractive,
      pointsx,pointsy,
      directionsx, directionsy,
      tangent_h,tangent_θ,
      nscans,nlos,logger
  )
end
