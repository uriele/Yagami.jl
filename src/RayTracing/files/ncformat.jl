using NCDatasets: NCDataset



"""
   $SIGNATURES

Create a `RayTracingProblem` from a Cairt file.
This struct represents a ray tracing problem with all the necessary data
to perform ray tracing calculations, including the atmosphere, refractive index, and satellite scan information.
# Arguments:
- `::Cairt`: A type representing the Cairt file format, which is used to read the ray tracing data.
- `filename::String`: The path to the Cairt file containing the ray tracing data.

# Key Arguments:
- `meantype::MT=GeometricMean()`: The mean type to be used for the refractive index calculation (default is `GeometricMean`).
- `model::AM=Ciddor()`: The air model to be used for the refractive index calculation (default is `Ciddor`).
- `earthmodel::EA=Fukushima()`: The earth approximation model to be used (default is `Fukushima`).
- `logger=nothing`: An optional logger to log information during the problem creation (default is `NullLogger`).
# Returns:
- `RayTracingProblem`: An instance of `RayTracingProblem` containing the data read from the Cairt file and the calculated refractive index and atmosphere.
"""
function NCRayTracingProblem(filename::String;
    meantype::MT=GeometricMean(), model::AM=Ciddor(),
    earthmodel::EA=Fukushima(),
    logger=nothing,
    humidity::Real=0.0, # used to override the humidity in the file
    co2ppm::Real=400.0, # used to override the CO2 concentration in the file
  )   where {MT<:MeanType,AM<:AirModel,
    EA<:EarthApproximation
  }
    nc=NCDataset(filename)

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
      @info "humidity: $(humidity*100)%, CO2 concentration: $co2ppm ppm"
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
      fn-> @. uconvert(Pa, fn*hPa)  |> # convert to Pa
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
      humidity=humidity,
      co2ppm=co2ppm,
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

    RayTracingProblem{T,MT,AM,EA,ATM,V,M,L}(
      filename, meantype, model, earthmodel,
      atmosphere,refractive,
      pointsx,pointsy,
      directionsx, directionsy,
      tangent_h,tangent_θ,
      nscans,nlos,logger )
end




Base.show(io::IO,::Type{<:RTP}) where {T<:AbstractFloat,MT<:MeanType,AM<:AirModel,
  EA<:EarthApproximation,ATM<:AtmosphereSetting,
  V<:AbstractVector{T},
  M<:AbstractMatrix{T},L<:AbstractLogger,RTP<:RayTracingProblem{T,MT,AM,EA,ATM,V,M,L}} = print(io, "RayTracingProblem{$T}")
Base.show(::IO, probl::RayTracingProblem{T,MT,AM,EA,ATM,V,M,L}
) where {T<:AbstractFloat,MT<:MeanType,AM<:AirModel,
  EA<:EarthApproximation,ATM<:AtmosphereSetting,
  V<:AbstractVector{T},
  M<:AbstractMatrix{T},L<:AbstractLogger} = begin
  println("RayTracingProblem{$T}:")
  println("   Filename: $(probl.filename)")
  println("   MeanType: $(MT)")
  println("   Air Model: $(AM)")
  println("   Earth Model: $(EA)")
  println("   Number of scans: $(probl.nscans)")
  println("   Number of lines of sight: $(probl.nlos)")
end
Base.propertynames(::RayTracingProblem) = PROBLPROPTOTAL

function Base.getproperty(probl::RayTracingProblem, sym::Symbol)
  if sym ∈ PROBLKNOTS
    return getfield(getfield(getfield(probl,:atmosphere),:temperature), sym)
  elseif sym ∈ PROBLATMPROP
    return getfield(getfield(probl,:atmosphere),sym)
  elseif sym ∈ PROBLPROPDIRECT
    return getfield(probl, sym)
  end
  throw(ArgumentError("Property $sym not found in RayTracingProblem."))
end

getlogger(probl::RayTracingProblem) = getfield(probl, :__logger)

############################################################################################
# helper functions
############################################################################################
@inline los_scan_to_index(prob::RayTracingProblem,los::Int,scan::Int) = begin
  nlos=prob.nlos
  nscan=prob.nscan
  @inline idx(los,scan)=if los<=nlos && scan<=nscan
    los+(scan-1)*(nlos)
  else
    error("Invalid los or scan index. los: $los, scan: $scan, nlos: $nlos, nscan: $nscan")
  end
  idx(los,scan)
end

EXPCAIRTHELPERS=[:wedge_refractive,:getlogger,
:num_wedges,:num_edges,:extrapolatetemperature,
:extrapolatepressure,:extrapolatehumidity,:extrapolateco2ppm,:extrapolatewavelength,
:getsatposition,:getsatdirection,:getquotetangent]

"""
   `wedge_refractive(problem::RayTracingProblem, i_theta::Int, j_h::Int;complement::Bool=false)`
Returns the refractive index at a given wedge edge specified by indices `i_theta` and `j_h` in a `RayTracingProblem`.
# Arguments:
- `problem::RayTracingProblem`: The ray tracing problem instance.
- `i_theta::Int`: The index of the theta knot (angle).
- `j_h::Int`: The index of the h knot (height).
# Key Arguments:
- `complement::Bool`: If true, returns the complement of the refractive index (n-1), that can be more useful in some cases (default is false).

# Returns:
- `Tuple{Tuple{Float64, Float64}, Float64}`: A tuple containing the theta and h values at the specified indices and the refractive index.
# Example:
```julia
julia> wedge_refractive(prob, 1, 1) # returns the refractive index for i_theta=1 and j_h=1
((0.0, 0.0), 1.0003)
```
"""
function wedge_refractive(problem::RayTracingProblem, i_theta::Int, j_h::Int;complement::Bool=false)
  knots_h = problem.temperature.knots_h
  knots_theta = problem.temperature.knots_θ[1:end-1]
  N, M = length(knots_h), length(knots_theta)
  if i_theta < 1 || i_theta > M || j_h < 1 || j_h > N
    error("Invalid indices: i_theta=$i_theta, j_h=$j_h, M=$M, N=$N")
  end

  n= problem.refractive[i_theta, j_h]
  if complement
    n = n - 1
  end
  return (knots_theta[i_theta], knots_h[j_h]), n
end

for value in (:temperature, :pressure, :humidity, :co2ppm, :wavelength)
  func= Symbol("node_", value)
  push!(EXPCAIRTHELPERS, func)
  string_doc = """
    `$func(prob::RayTracingProblem, i_theta::Int, j_h::Int)`
    Returns the nodes of the $value for a given `RayTracingProblem` at specified wedge edge.
    # Arguments:
    - `prob::RayTracingProblem`: The ray tracing problem instance.
    - `i_theta::Int`: The index of the theta knot (angle).
    - `j_h::Int`: The index of the h knot (height).
    # Returns:
    - `Tuple{Int,Int},T`: A tuple containing the theta and h values at the specified indices and the value of $value.

    # Example:
    ```julia
    julia> $func(prob, 1, 1) # returns the nodes for i_theta=1 and j_h=1
    (0.0, 0.0, 273.15)
    ```
    """
    @eval begin
    @doc $string_doc
    @inline function $func(problem::RayTracingProblem, i_theta::Int, j_h::Int)
      knots_h = problem.temperature.knots_h
      knots_theta = problem.temperature.knots_θ[1:end-1]
      N, M = length(knots_h), length(knots_theta)
      if i_theta < 1 || i_theta > M || j_h < 1 || j_h > N
        error("Invalid indices: i_theta=$i_theta, j_h=$j_h, M=$M, N=$N")
      end
      return (knots_theta[i_theta], knots_h[j_h]), problem.$value(float(knots_theta[i_theta]), float(knots_h[j_h]))
    end
  end
end

"""
   `numwedges(prob::RayTracingProblem)`
Returns the number of wedges in a `RayTracingProblem` as N x M, where N is the number of angle knots and M is the number of height knots.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
# Returns:
- `Tuple{Int, Int}`: A tuple containing the number of wedges (N, M).

## Note:
due to the circularity of the matrix along the theta dimension, the number of wedges in the theta dimension is the same
as the number of edges in the theta direction, while the number of wedges in the height dimension is one less than the number of edges in the h direction.

# Example:
```julia
julia> numwedges(prob) # returns the number of wedges in the problem
(10, 20)
```
See also: [`num_edges`](@ref)
"""
num_wedges(prob::RayTracingProblem) = size(prob.refractive,1)-1,size(prob.refractive,2)
"""
   `num_edges(prob::RayTracingProblem)`
Returns the number of edges in a `RayTracingProblem` as N x M, where N is the number of angle knots and M is the number of height knots.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
# Returns:
- `Tuple{Int, Int}`: A tuple containing the number of edges (N, M).
# Note:
Due to the circularity of the matrix along the theta dimension, the number of edges in the theta direction is one less than the number of knots in the theta direction, while the number of edges in the height dimension is one less than the number of knots in the height direction.
# Example:
```julia
julia> num_edges(prob) # returns the number of edges in the problem
(10, 19)
```
See also: [`num_wedges`](@ref)
"""
num_edges(prob::RayTracingProblem) = size(prob.temperature.knots_θ,1)-1,size(prob.temperature.knots_h,1)


"""
   `extrapolatetemperature(prob::RayTracingProblem, θ::AbstractFloat, h::AbstractFloat)`
   `extrapolatetemperature(prob::RayTracingProblem; theta=0.0, h=0.0)`
Extrapolates the temperature at a given line of angle `θ` and height `h` for a given `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `θ::AbstractFloat`: The angle in degree and geodesic coordinates.
- `h::AbstractFloat`: The altitude with respect to the surface in km.
# Returns:
- `Float64`: The extrapolated temperature in Kelvin

# Key Arguments:
- `theta::AbstractFloat`: The angle in degrees (default is 0.0).
- `h::AbstractFloat`: The height in kilometers (default is 0.0).

# Example:
```julia
julia> extrapolatetemperature(prob, 10.0, 5.0) # angle in degrees, height in km
273.15
julia> extrapolatetemperature(prob; theta=10.0, h=5.0) # using keyword arguments
275.15
julia> extrapolatetemperature(prob) # defaults to theta=0.0, h=0.0
278.15
```
"""
extrapolatetemperature(prob::RTP,θ::T,h::T) where {T<:AbstractFloat,RTP<:RayTracingProblem{T}}=prob.temperature(float(θ),float(h))
function extrapolatetemperature(prob::RTP;theta::Real=0.0,h::Real=0.0) where {RTP<:RayTracingProblem}
  T=eltype(theta)
  extrapolatetemperature(T(prob),T(theta),T(h))
end
"""
   `extrapolatepressure(prob::RayTracingProblem, θ::AbstractFloat, h::AbstractFloat)`
    `extrapolatepressure(prob::RayTracingProblem; theta=0.0, h=0.0)`
Extrapolates the pressure at a given line of angle `θ` and height `h` for a given `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `θ::AbstractFloat`: The angle in degree and geodesic coordinates.
- `h::AbstractFloat`: The altitude with respect to the surface in km.
# Returns:
- `Float64`: The extrapolated pressure in Pa

# Key Arguments:
- `theta::AbstractFloat`: The angle in degrees (default is 0.0).
- `h::AbstractFloat`: The height in kilometers (default is 0.0).

# Example:
```julia
julia
> extrapolatepressure(prob, 10.0, 5.0) # angle in degrees, height in km
101325.0
julia> extrapolatepressure(prob; theta=10.0, h=5.0) # using keyword arguments
101325.0
julia> extrapolatepressure(prob) # defaults to theta=0.0, h=0.0
101325.0
```

"""
extrapolatepressure(prob::RTP,θ::T,h::T) where {T<:AbstractFloat,RTP<:RayTracingProblem{T}} =prob.pressure(float(θ),float(h))
function extrapolatepressure(prob::RTP;theta::Real=0.0,h::Real=0.0) where {RTP<:RayTracingProblem}
  T=eltype(theta)
  extrapolatepressure(T(prob),T(theta),T(h))
end

"""
    `extrapolatehumidity(prob::RayTracingProblem, θ::AbstractFloat, h::AbstractFloat)`
    `extrapolatehumidity(prob::RayTracingProblem; theta=0.0, h=0.0)`
Extrapolates the humidity at a given line of angle `θ` and height `h` for a given `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `θ::AbstractFloat`: The angle in degree and geodesic coordinates.
- `h::AbstractFloat`: The altitude with respect to the surface in km.
# Returns:
- `Float64`: The extrapolated humidity in %
"""
extrapolatehumidity(prob::RTP,θ::T,h::T) where {T<:AbstractFloat,RTP<:RayTracingProblem{T}} =prob.humidity(float(θ),float(h))*100
function extrapolatehumidity(prob::RTP;theta::Real=0.0,h::Real=0.0) where {RTP<:RayTracingProblem}
  T=eltype(theta)
  extrapolatehumidity(T(prob),T(theta),T(h))
end


"""
    `extrapolateco2ppm(prob::RayTracingProblem, θ::AbstractFloat, h::AbstractFloat)`
    `extrapolateco2ppm(prob::RayTracingProblem; theta=0.0, h=0.0)`
Extrapolates the CO2 concentration in ppm at a given line of angle `θ` and height `h` for a given `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `θ::AbstractFloat`: The angle in degree and geodesic coordinates.
- `h::AbstractFloat`: The altitude with respect to the surface in km.
# Returns:
- `Float64`: The extrapolated CO2 concentration in ppm

# Key Arguments:
- `theta::AbstractFloat`: The angle in degrees (default is 0.0).
- `h::AbstractFloat`: The height in kilometers (default is 0.0).

# Example:
```julia
julia> extrapolateco2ppm(prob, 10.0, 5.0) # angle in degrees, height in km
400.0
julia> extrapolateco2ppm(prob; theta=10.0, h=5.0) # using keyword arguments
400.0
julia> extrapolateco2ppm(prob) # defaults to theta=0.0, h=0.0
400.0
```
"""
extrapolateco2ppm(prob::RTP,θ::T,h::T) where {T<:AbstractFloat,RTP<:RayTracingProblem{T}} =prob.co2ppm(float(θ),float(h))
function extrapolateco2ppm(prob::RTP;theta::Real=0.0,h::Real=0.0) where {RTP<:RayTracingProblem}
  T=eltype(theta)
  extrapolateco2ppm(T(prob),T(theta),T(h))
end
"""
    `extrapolatewavelength(prob::RayTracingProblem, θ::AbstractFloat, h::AbstractFloat)`
    `extrapolatewavelength(prob::RayTracingProblem; theta=0.0, h=0.0)`
Extrapolates the wavelength at a given line of angle `θ` and height `h` for a given `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `θ::AbstractFloat`: The angle in degree and geodesic coordinates.
- `h::AbstractFloat`: The altitude with respect to the surface in km.
# Returns:
- `Float64`: The extrapolated wavelength in μm

# Key Arguments:
- `theta::AbstractFloat`: The angle in degrees (default is 0.0).
- `h::AbstractFloat`: The height in kilometers (default is 0.0).

# Example:
```julia
julia> extrapolatewavelength(prob, 10.0, 5.0) # angle in degrees, height in km
1.0
julia> extrapolatewavelength(prob; theta=10.0, h=5.0) # using keyword arguments
1.0
julia> extrapolatewavelength(prob) # defaults to theta=0.0, h=0.0
1.0
```
"""
extrapolatewavelength(prob::RTP,θ::T,h::T) where {T<:AbstractFloat,RTP<:RayTracingProblem{T}} =prob.wavelength(float(θ),float(h))
extrapolatewavelength(prob::RTP;theta::Real=0.0,h::Real=0.0) where {T<:AbstractFloat,RTP<:RayTracingProblem{T}}=extrapolatewavelength(T(prob),T(theta),T(h))

"""
  `getsatposition(prob::RayTracingProblem, los::Int, scan::Int)`
  `getsatposition(prob::RayTracingProblem; los=1, scan=1)`
Returns the satellite initial position for a given line of sight (los) and scan index in a `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `los::Int`: The line of sight index.
- `scan::Int`: The scan index.
# Returns:
- `Vector{Float64}`: A vector containing the x and y coordinates of the satellite position.

# Key Arguments:
- `los::Int`: The line of sight index (default is 1).
- `scan::Int`: The scan index (default is 1).

# Example:
```julia
julia
> getsatposition(prob, 1, 1) # returns the position for los=1 and scan=1
[100.0, 200.0]
julia> getsatposition(prob; los=1, scan=1) # using keyword arguments
[100.0, 200.0]
julia> getsatposition(prob) # defaults to los=1, scan=1
[100.0, 200.0]
"""
getsatposition(prob::RayTracingProblem,los::Int,scan::Int) = begin
  ii=los_scan_to_index(prob,los,scan)
  return [prob.pointx[ii], prob.pointy[ii]]
end
getsatposition(prob::RayTracingProblem;los::Int=1,scan::Int=1) = getsatposition(prob,los,scan)


"""
  `getsatdirection(prob::RayTracingProblem, los::Int, scan::Int)`
  `getsatdirection(prob::RayTracingProblem; los=1, scan=1)`
Returns the satellite direction (from the nadir angle) for a given line of sight (los) and scan index in a `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `los::Int`: The line of sight index.
- `scan::Int`: The scan index.
# Returns:
- `Vector{Float64}`: A vector containing the x and y components of the satellite direction.

# Key Arguments:
- `los::Int`: The line of sight index (default is 1).
- `scan::Int`: The scan index (default is 1).
# Example:
```julia
julia> getsatdirection(prob, 1, 1) # returns the direction for los=1 and scan=1
[0.1, 0.2]
julia> getsatdirection(prob; los=1, scan=1) # using keyword arguments
[0.1, 0.2]
julia> getsatdirection(prob) # defaults to los=1, scan=1
[0.1, 0.2]
````
"""
getsatdirection(prob::RayTracingProblem,los::Int,scan::Int) = begin
  ii=los_scan_to_index(prob,los,scan)
  return [prob.directionx[ii], prob.directiony[ii]]
end
getsatdirection(prob::RayTracingProblem;los::Int=1,scan::Int=1) = getsatdirection(prob,los,scan)

"""
  `getquotetangent(prob::RayTracingProblem, los::Int, scan::Int)`
  `getquotetangent(prob::RayTracingProblem; los=1, scan=1)`
Returns the tangent point coordinates for a given line of sight (los) and scan index in a `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `los::Int`: The line of sight index.
- `scan::Int`: The scan index.
# Returns:
- `Tuple{Float64, Float64}`: A tuple containing the tangent height and azimuth angle at the specified los and scan.
# Key Arguments:
- `los::Int`: The line of sight index (default is 1).
- `scan::Int`: The scan index (default is 1).
# Example:
```julia
julia> getquotetangent(prob, 1, 1) # returns the tangent point for los=1 and scan=1
(100.0, 45.0)
julia> getquotetangent(prob; los=1, scan=1) # using keyword arguments
(100.0, 45.0)
julia> getquotetangent(prob) # defaults to los=1, scan=1
(100.0, 45.0)
```
"""
getquotetangent(prob::RayTracingProblem,los::Int,scan::Int) = begin
  ii=los_scan_to_index(prob,los,scan)
  return (prob.tangent_h[ii], prob.tangent_θ[ii])
end
getquotetangent(prob::RayTracingProblem;los::Int=1,scan::Int=1) = getquotetangent(prob,los,scan)
