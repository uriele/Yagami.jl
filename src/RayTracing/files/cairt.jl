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



struct RayTracingProblem{T,MT<:MeanType,AM<:AirModel,EA<:EarthApproximation,ATM<:AtmosphereSetting,V<:AbstractVector{T},M<:AbstractMatrix{T},L<:AbstractLogger}
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
  __logger::L


  function RayTracingProblem(filename::String;
    meantype::MT=GeometricMean(), model::AM=Ciddor(),
    earthmodel::EA=Fukushima(),
    logger=nothing,
    ) where {MT<:MeanType,AM<:AirModel,
    EA<:EarthApproximation}
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
   `extrapolatetemperature(prob::RayTracingProblem, θ::Real, h::Real)`
   `extrapolatetemperature(prob::RayTracingProblem; theta=0.0, h=0.0)`
Extrapolates the temperature at a given line of angle `θ` and height `h` for a given `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `θ::Real`: The angle in degree and geodesic coordinates.
- `h::Real`: The altitude with respect to the surface in km.
# Returns:
- `Float64`: The extrapolated temperature in Kelvin

# Key Arguments:
- `theta::Real`: The angle in degrees (default is 0.0).
- `h::Real`: The height in kilometers (default is 0.0).

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
extrapolatetemperature(prob::RayTracingProblem,θ::Real,h::Real)=prob.temperature(float(θ),float(h))
extrapolatetemperature(prob::RayTracingProblem;theta::Real=0.0,h::Real=0.0)=extrapolatetemperature(prob,theta,h)
"""
   `extrapolatepressure(prob::RayTracingProblem, θ::Real, h::Real)`
    `extrapolatepressure(prob::RayTracingProblem; theta=0.0, h=0.0)`
Extrapolates the pressure at a given line of angle `θ` and height `h` for a given `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `θ::Real`: The angle in degree and geodesic coordinates.
- `h::Real`: The altitude with respect to the surface in km.
# Returns:
- `Float64`: The extrapolated pressure in Pa

# Key Arguments:
- `theta::Real`: The angle in degrees (default is 0.0).
- `h::Real`: The height in kilometers (default is 0.0).

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
extrapolatepressure(prob::RayTracingProblem,θ::Real,h::Real)=prob.pressure(float(θ),float(h))
extrapolatepressure(prob::RayTracingProblem;theta::Real=0.0,h::Real=0.0)=extrapolatepressure(prob,theta,h)

"""
    `extrapolatehumidity(prob::RayTracingProblem, θ::Real, h::Real)`
    `extrapolatehumidity(prob::RayTracingProblem; theta=0.0, h=0.0)`
Extrapolates the humidity at a given line of angle `θ` and height `h` for a given `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `θ::Real`: The angle in degree and geodesic coordinates.
- `h::Real`: The altitude with respect to the surface in km.
# Returns:
- `Float64`: The extrapolated humidity in %
"""
extrapolatehumidity(prob::RayTracingProblem,θ::Real,h::Real)=prob.humidity(float(θ),float(h))*100
extrapolatehumidity(prob::RayTracingProblem;theta::Real=0.0,h::Real=0.0)=extrapolatehumidity(prob,theta,h)


"""
    `extrapolateco2ppm(prob::RayTracingProblem, θ::Real, h::Real)`
    `extrapolateco2ppm(prob::RayTracingProblem; theta=0.0, h=0.0)`
Extrapolates the CO2 concentration in ppm at a given line of angle `θ` and height `h` for a given `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `θ::Real`: The angle in degree and geodesic coordinates.
- `h::Real`: The altitude with respect to the surface in km.
# Returns:
- `Float64`: The extrapolated CO2 concentration in ppm

# Key Arguments:
- `theta::Real`: The angle in degrees (default is 0.0).
- `h::Real`: The height in kilometers (default is 0.0).

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
extrapolateco2ppm(prob::RayTracingProblem,θ::Real,h::Real)=prob.co2ppm(float(θ),float(h))
extrapolateco2ppm(prob::RayTracingProblem;theta::Real=0.0,h::Real=0.0)=extrapolateco2ppm(prob,theta,h)
"""
    `extrapolatewavelength(prob::RayTracingProblem, θ::Real, h::Real)`
    `extrapolatewavelength(prob::RayTracingProblem; theta=0.0, h=0.0)`
Extrapolates the wavelength at a given line of angle `θ` and height `h` for a given `RayTracingProblem`.
# Arguments:
- `prob::RayTracingProblem`: The ray tracing problem instance.
- `θ::Real`: The angle in degree and geodesic coordinates.
- `h::Real`: The altitude with respect to the surface in km.
# Returns:
- `Float64`: The extrapolated wavelength in μm

# Key Arguments:
- `theta::Real`: The angle in degrees (default is 0.0).
- `h::Real`: The height in kilometers (default is 0.0).

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
extrapolatewavelength(prob::RayTracingProblem,θ::Real,h::Real)=prob.wavelength(float(θ),float(h))
extrapolatewavelength(prob::RayTracingProblem;theta::Real=0.0,h::Real=0.0)=extrapolatewavelength(prob,theta,h)

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
