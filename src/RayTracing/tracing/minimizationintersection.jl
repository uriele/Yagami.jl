const ACCEPTABLE_TOLERANCE::Float64 =1e-5
const FREESPACE::Float64 = 1.0 # free space index
struct NoFunc end
abstract type AbstractResult{T<:AbstractFloat} end
struct SimpleResult{T} <:AbstractResult{T}
  pointx::T
  pointy::T
  directionx::T
  directiony::T
  altitude::T
  azimuth::T
  length::T
  i::Int
  j::Int
end
SimpleResult(pointx::T=MAJORAXIS(), pointy::T=MINORAXIS(), directionx::T=0.0, directiony::T=-1.0,
  altitude::T=T(0), azimuth::T=T(0), length::T=T(0), i::Int=0, j::Int=0) where {T<:AbstractFloat} =
  SimpleResult{T}(pointx, pointy, directionx, directiony, altitude, azimuth, length, i, j)

SimpleResult()
# pass traits of the Earth Approximation
EarthApproximation(x) = EarthApproximation(typeof(x)) # convert to type
EarthApproximation(::Type)= NoFunc() # no type conversion
EarthApproximation(::T) where {T<:Type{<:Fukushima}} = Fukushima
EarthApproximation(::T) where {T<:Type{<:Bowring}} = Bowring

@inline earthmodelfunction(x)=earthmodelfunction(EarthApproximation(x))
@inline earthmodelfunction(::NoFunc) = error("No Earth Approximation specified. Use Fukushima or Bowring.")
@inline earthmodelfunction(::Type{<:Fukushima}) = ___altitudeangle_fukushima

@inline earthmodelfunction(::Type{<:Bowring}) = ___altitudeangle_bowring

@inline function __update_ray!(zb::Zbrent,pointx::T,pointy::T,directionx::T,directiony::T,
  hmin::T,hmax::T,θmin::T,θmax::T,i::Int,j::Int,isdescending::Bool,islevel::Bool) where {T<:AbstractFloat}
  __setpoint!(zb, pointx, pointy) # set the new point position
  __setdirection!(zb, directionx, directiony) # set the new direction
  __sethlims!(zb, hmin, hmax) # set the new h limits
  __setθlims!(zb, θmin, θmax) # set the new θ limits
  __seti!(zb, i) # set the new i index
  __setj!(zb, j) # set the new j index
  __setislevel!(zb, islevel) # set to islevel
  __setdescending!(zb, isdescending) # reset to descending
  __setbracket!(zb, 0.0, 0.0, 0.0) # reset the bracket
end

@inline function __update_results!(res::A,iter::Int, pointx::T,pointy::T,
  directionx::T,directiony::T,
  altitude::T,azimuth::T,
  len_t::T,i::Int,j::Int) where {T<:AbstractFloat,R<:AbstractResult{T},A<:AbstractVector{R}}

  res.pointx[iter] = pointx
  res.pointy[iter] = pointy
  res.directionx[iter] = directionx
  res.directiony[iter] = directiony
  res.altitude[iter] = altitude
  res.azimuth[iter] = azimuth
  res.length[iter] = len_t
  res.i[iter] = i
  res.j[iter] = j
end

@inline function logstep(logger,
  iter::Int, rx::T, ry::T,
  dirx::T, diry::T, len_t::T, idx_i::Int, idx_j::Int,
  hk::T, θc::T,isdecending::Bool) where {T<:AbstractFloat}
  with_logger(logger) do
    @info "$(numshort(iter))  $(numshort(idx_i))  $(numshort(idx_j))  $(numshort(rx))  $(numshort(ry))  $(numshort(dirx))  $(numshort(diry))  $(numshort(len_t))  $(numshort(hk))  $(numshort(θc))  $(isdecending ? "descending" : "ascending")"
  end
end

@inline __innerfunc(z::Zbrent, x::T,y::T) where T = getfield(z.f, :__f)(x,y,MAJORAXIS(T),MINORAXIS(T),COMPLECCENTRICITY²(T))

function raytracing!(
  results::RR,
  problem::RP,
  itermax::Int = 100,
  tol::T = T(1e-10)
) where {
  T<:AbstractFloat,
  RE<:AbstractResult{T},
  RP<:RayTracingProblem,
  RR<:AbstractMatrix{RE}
}
  # Unpack the problem
  pointsx = problem.pointsx
  pointsy = problem.pointsy
  directionsx = problem.directionsx
  directionsy = problem.directionsy
  knots_h = problem.atmosphere.temperature.knots_h
  knots_θ = problem.atmosphere.temperature.knots_θ
  refractive_grid = problem.refractive

  model = problem.earthmodel

  logger = getlogger(problem)
  nlos = problem.nlos
  nscans = problem.nscans
  with_logger(logger) do
    @info "============================================================"
    @info "Ray Tracing Problem"
    @info "============================================================"
    @info "    $(textshort("iter")) $(textshort("i")) $(textshort("j")) $(textshort("x [km]")) $(textshort("y [km]")) $(textshort("dX")) $(textshort("dY")) $(textshort("h [km]")) $(textshort("θ [°]"))"
  end


  zb = Zbrent(T,DistanceFunc(earthmodelfunction(model)),itermax,tol)
  @inbounds for i in eachindex(pointsx)

    with_logger(logger) do
      scan= div(i-1,nscans)+ 1
      los= mod1(i,nlos)
      @info "scan: $(numshort(scan))  los: $(numshort(los)) "
    end
    @inline __solve!(view(results,:,i),zb, pointsx[i], pointsy[i], directionsx[i], directionsy[i],view(knots_h,:), view(knots_θ,:), view(refractive_grid,:,:),logger)
  end
end



@inline function __solve!(res::A, zb::Zbrent{F,T}, pointx::T,pointy::T,directionx::T,directiony::T,
  knots_h::V,knots_θ::V,refractive_grid::MA,logger::L=NullLogger()) where {F,T<:AbstractFloat,R<:AbstractResult{T},A<:AbstractVector{R},V<:AbstractVector{T},MA<:AbstractMatrix{T},L<:AbstractLogger}

  N= length(knots_θ)-1
  M= length(knots_h)-1
  niters= length(res)
  # Initialize the ray structure
  __setpoint!(zb, pointx, pointy)
  __setdirection!(zb, directionx, directiony)
  hmax,_= __getinnerf(zb, pointx, pointy)
  hmin = knots_h[1]

  θmax = Inf
  θmin = -Inf

  n_i = FREESPACE

  __sethlims!(zb, hmin, hmax) # set the height limits
  __setθlims!(zb, -10., 10.)
  ###########################
  __setdescending!(zb) # reset to descending
  @inbounds for iter in 1:niters
    # find the minimium in the bracket

    # for caching
    rx = getpointx(zb)
    ry = getpointy(zb)
    dirx = getdirectionx(zb)
    diry = getdirectiony(zb)
    idx_i= geti(zb)
    idx_j= getj(zb)
    islevel = true # is the ray level? assume true every iteration

    isdescending = getdescending(zb)

    bracketmin(zb)

    Δ = abs(zb.b - zb.a) # get the difference between the bounds
    # update the bracket if minimum found
    if Δ< 10*zb.tol || zb.b<zb.a
      @info "iter: $(iter) - Bracket too small, resetting to descending"
      isdescending=false # reset to descending
      bracketmin(zb) # recalculate the bracket
    end
    # solve
    findmin(zb) # find the minimum

    len_t = zb.x # get the length of the ray
    fmin  = zb.fx
    # trial
    _x = rx+len_t * dirx
    _y = ry+len_t * diry
    hk,θk = __getinnerf(zb, _x, _y) # get the height and angle at the point
    # check if it is in the range (θmin, θmax)
    θc = clampangle(θk, θmin, θmax)

    # normal direction with respect to the intersection point
    # outward normal for the descending ray
    normalx = cosd(θc)
    normaly = sind(θc)

    if iter == 1
      if !isdescending
         break # if not descending, break the loop it means that it never entered the atmosphere
      end
      idx_j = 1
      idx_i = searchsortedlast(view(knots_θ,1:N), θc)
      # we know the ray is descending already
    elseif θc == θmin || θc == θmax
      # Find the correct value of t by using the intersection between two lines
      islevel = false
      # Find the correct value of t by using the intersection between two lines
      N₀ = MAJORAXIS(T)/sqrt(1-ECCENTRICITY²(T)*normaly*normaly) # normal to the Earth surface
      earthx = N₀*normalx
      earthy = COMPLECCENTRICITY²(T)*N₀*normaly
      len_t,hk = intersectionrayray(rx, ry, dirx, diry, earthx, earthy, normalx, normaly)
      # important: this is the direction for θc=θmax
      normalx,normaly = normaly, -normalx # swap the normal direction
      isθmax = θc == θmax

      if isθmax
        idx_i = mod1(idx_i+1,N) # set the next index
      else
        idx_i = mod1(idx_i-1,N) # set the previous index
        normalx,normaly = -normalx, -normaly # reverse the normal direction
      end
    elseif abs(fmin)< ACCEPTABLE_TOLERANCE
      if getdescending(zb)
        idx_j += 1 # set the next index
      else
        idx_j -= 1 # set the previous index
        normalx,normaly = -normalx, -normaly # reverse the normal direction
      end
    end

    # check if still in the atmosphere otherwise break the loop
    if (idx_j < 1 || idx_j >= M)
      break # if the index is out of bounds, break the loop
    end
    # update the point position
    rx +=  len_t * dirx
    ry +=  len_t * diry
    ########################################################################################

    n_transmitted = refractive_grid[idx_i,idx_j] # get the refractive index at the current
    # update the direction
    dirx,diry = snellslaw(normalx, normaly, dirx, diry, n_i, n_transmitted)

    # update the refractive index for the next iteration
    n_i = n_transmitted # update the refractive index for the next iteration

    # h in descending order
    hmax = knots_h[idx_j]
    hmin = knots_h[idx_j+1] # get the hmax and hmin from the knots
    # note: knots_θ has an extra element to account for the last angle
    # theta in ascending order
    θmin = knots_θ[idx_i]
    θmax = knots_θ[idx_i+1] # knots_θ has an extra element to account for the last angle
    # update for next iteration

    __update_ray!(zb, rx, ry, dirx, diry,
      hmin, hmax, θmin, θmax,
      idx_i, idx_j,
      isdescending, islevel) # update the ray

    # update the results
    #######################################################################################
    __update_results!(res, iter, rx, ry, dirx, diry,
      hk, θc, len_t, idx_i, idx_j) # update the results
    ########################################################################################
    # log the step
    with_logger(logger) do
      logstep(logger, iter, rx, ry, dirx, diry, len_t, idx_i, idx_j, hk, θc,isdescending)
    end
  end
end
