
# pass traits of the Earth Approximation
EarthApproximation(x) = EarthApproximation(typeof(x)) # convert to type
EarthApproximation(::Type)= NoFunc() # no type conversion
EarthApproximation(::T) where {T<:Type{<:Fukushima}} = Fukushima
EarthApproximation(::T) where {T<:Type{<:Bowring}} = Bowring

@inline earthmodelfunction(x)=earthmodelfunction(EarthApproximation(x))
@inline earthmodelfunction(::NoFunc) = error("No Earth Approximation specified. Use Fukushima or Bowring.")
@inline earthmodelfunction(::Type{<:Fukushima}) = ___altitudeangle_fukushima

@inline earthmodelfunction(::Type{<:Bowring}) = ___altitudeangle_bowring

@inline function __update_ray!(zb::Z,len_t::T,pointx::T,pointy::T,
  directionx::T,directiony::T,
  hmin::T,hmax::T,θmin::T,θmax::T,
  i::Int,j::Int,n_i::T,isdescending::Bool,islevel::Bool
) where {F,T<:AbstractFloat,Z<:Zbrent{F,T}}
  __setpoint!(zb, pointx, pointy) # set the new point position
  __setdirection!(zb, directionx, directiony) # set the new direction
  __sethlims!(zb, hmin, hmax) # set the new h limits
  __setθlims!(zb, θmin, θmax) # set the new θ limits
  __seti!(zb, i) # set the new i index
  __setj!(zb, j) # set the new j index
  __setn!(zb, n_i) # set the new refractive index
  __setislevel!(zb, islevel) # set to islevel
  __setdescending!(zb, isdescending) # reset to descending
  __setbracket!(zb, len_t, len_t, len_t) # reset the bracket
end



@inline function __update_results!(res::A,iter::Int,
  zb::Z,n_i::T=FREESPACE
) where {F,T<:AbstractFloat,R<:AbstractResult{T},
  A<:AbstractVector{R},Z<:Zbrent{F,T}
}
  pointx,pointy = getpointx(zb), getpointy(zb)
  res.pointx[iter+1] = pointx
  res.pointy[iter+1] = pointy
  res.directionx[iter+1] = getdirectionx(zb)
  res.directiony[iter+1] = getdirectiony(zb)
  if iter>0
    res.length_t[iter] = zb.x
  end
  res.i[iter+1] = geti(zb)
  res.j[iter+1] = getj(zb)
  altitude, azimuth = __getinnerf(zb, pointx,pointy) # get the altitude and azimuth at the point
  # More efficient ways exist but requires to update the Zbrent structure
  res.altitude[iter+1] = altitude
  res.azimuth[iter+1] = azimuth
  res.islevel[iter+1] = getislevel(zb) # is the ray level?
  res.descending[iter+1] = getdescending(zb) # is the ray descending?

end

@inline __innerfunc(z::Z, x::T,y::T) where {F,T,Z<:Zbrent{F,T}} = getfield(z.f, :__f)(x,y,MAJORAXIS(T),MINORAXIS(T),COMPLECCENTRICITY²(T))

function raytracing!(
  results::RR,
  problem::RP,
  itermax::Int = 100,
  tol::T = convert(T,1e-10)
) where {
  T<:AbstractFloat,
  RE<:AbstractResult{T},
  RP<:RayTracingProblem{T},
  RR<:AbstractMatrix{RE}
}
  # Unpack the problem
  pointsx = problem.pointsx
  pointsy = problem.pointsy
  directionsx = problem.directionsx
  directionsy = problem.directionsy
  refractive_grid = view(problem.refractive,:,:)
  Ni,Mi= size(refractive_grid) # number of angles and heights
  knots_h = view(problem.atmosphere.temperature.knots_h,1:Ni)
  knots_θ = view(problem.atmosphere.temperature.knots_θ,:)
  refractive_grid = problem.refractive

  model = problem.earthmodel

  logger = getlogger(problem)
  nlos = problem.nlos
  nscans = problem.nscans

  zb = Zbrent(T,DistanceFunc(earthmodelfunction(model)),itermax,tol)
  @inbounds for i in eachindex(pointsx)
    @inline __solve!(view(results,:,i),zb, pointsx[i], pointsy[i], directionsx[i], directionsy[i],knots_h, knots_θ, refractive_grid)
  end

  write_tracing_log(logger, results,refractive_grid)
end



function raytracing_parallel!(
  results::RR,
  problem::RP,
  itermax::Int = 100,
  tol::T = convert(T,1e-10)
) where {
  T<:AbstractFloat,
  RE<:AbstractResult{T},
  RP<:RayTracingProblem{T},
  RR<:AbstractMatrix{RE}
}
  # Unpack the problem
  pointsx = problem.pointsx
  pointsy = problem.pointsy
  directionsx = problem.directionsx
  directionsy = problem.directionsy
  refractive_grid = view(problem.refractive,:,:)
  Ni,Mi= size(refractive_grid) # number of angles and heights
  knots_h = view(problem.atmosphere.temperature.knots_h,1:Ni)
  knots_θ = view(problem.atmosphere.temperature.knots_θ,:)
  refractive_grid = problem.refractive

  model = problem.earthmodel

  logger = getlogger(problem)
  nlos = problem.nlos
  nscans = problem.nscans

  zb = Zbrent(T,DistanceFunc(earthmodelfunction(model)),itermax,tol)

  nthreads=Threads.nthreads()
  zb_buffer = Vector{typeof(zb)}(undef, nthreads)
  @inbounds for i in eachindex(zb_buffer)
    zb_buffer[i] = deepcopy(zb) # create a copy of the zbrent structure for each thread
  end

  @batch for i in eachindex(pointsx)
    @inbounds @inline __solve!(view(results,:,i),zb_buffer[Threads.threadid()], pointsx[i], pointsy[i], directionsx[i], directionsy[i],knots_h, knots_θ, refractive_grid)
  end

  write_tracing_log(logger, results,refractive_grid)
end



#initialize the solver
@inline function solveinit!(zb::ZB,n_i::T,pointx::T,pointy::T,
  directionx::T,directiony::T,hmin::T=zero(T)
) where {F,T<:AbstractFloat,ZB<:Zbrent{F,T}}
  # Initialize the Zbrent structure
  hmax,_= __getinnerf(zb, pointx, pointy)
  __setpoint!(zb, pointx, pointy)
  __setdirection!(zb, directionx, directiony)
  __sethlims!(zb, hmin, hmax) # set the height limits to zero
  __setθlims!(zb, T(-Inf), T(Inf)) # set the angle limits
  __seti!(zb, 0) # set the i index
  __setj!(zb, 0) # set the j index
  __setn!(zb, n_i) # set the refractive index
  __setislevel!(zb, true) # set to islevel
  __setdescending!(zb) # reset to descending
  __setbracket!(zb, 0.0, 0.0, 0.0) # reset the bracket
end

#inplace solve function
function solvenext!(iter::Int,zb::Z,
  knots_h::Vh,knots_θ::Vθ,refractive_grid::MA
)::Bool where {F,T<:AbstractFloat,Z<:Zbrent{F,T},
  Vh<:AbstractVector{T},Vθ<:AbstractVector{T},
  MA<:AbstractMatrix{T}
}
    # caching the values
    ###########################################################
    # earth model parameters
    ###########################################################
    a =MAJORAXIS(T) # major axis
    a²=MAJORAXIS(T)*MAJORAXIS(T) # major axis squared
    b²=MINORAXIS(T)*MINORAXIS(T) # minor axis squared
    c²=COMPLECCENTRICITY²(T) # complement of the eccentricity squared
    e² = ECCENTRICITY²(T) # eccentricity squared
    ###########################################################
    tol= zb.tol # tolerance
    rx = getpointx(zb)
    ry = getpointy(zb)
    dirx = getdirectionx(zb)
    diry = getdirectiony(zb)
    n_i  = getn(zb) # refractive index at the current point
    idx_i= geti(zb)
    idx_j= getj(zb)
    islevel = true # is the ray level? assume true every iteration
    isdescending = getdescending(zb)
    θmin = getθmin(zb) # minimum angle
    θmax = getθmax(zb) # maximum angle
    hmin = gethmin(zb) # minimum height
    hmax = gethmax(zb) # maximum height
    N,M  = size(refractive_grid) # number of angles and heights

    bracketmin(zb,zero(T))

    # Try to find good condition for switching
    flag1= zb.b <zb.a
    #flag2= abs(zb.a-zb.b)< tol
    # update the bracket if minimum found
    if flag1
      isdescending=!isdescending # reset to descending
      __setdescending!(zb, isdescending) # update the descending flag
      # try to use the angle to kick it out of the minimum
      h_curr,θ_curr = __getinnerf(zb, rx+tol*dirx, ry+tol*diry) # get the height and angle at the point

      dcurr = isdescending ? h_curr - hmin : hmax - h_curr # distance to the current height

      # check if the distance is too small
      dcurr < 0 && return false # either the levels are too close or there is some bigger issue

      # create an adhoc bracket to escape the maximum
      dminθ = mod(θ_curr - θmin,360) # distance to the minimum angle
      dmaxθ = mod(θmax - θ_curr,360) # distance to the maximum angle
      θkick = θ_curr+max(dminθ, dmaxθ) # choose the angle with the maximum distance
      cosθkick = cosd(θkick) # kick direction in x
      sinθkick = sind(θkick) # kick direction in y
      N₀= a/sqrt(1-e²*sinθkick*sinθkick) # normal to the Earth surface
      _pointx = N₀ * cosθkick # kick point in x
      _pointy = c² * N₀ * sinθkick # kick point in y
      _dirx   = _pointx/a²
      _diry   = _pointy/b² # kick direction in y
      _normdir = hypot(_dirx, _diry) # normalize the kick direction
      _dirx /= _normdir
      _diry /= _normdir # normalize the kick direction
       x, _ = intersectionrayray(rx, ry, dirx, diry, _pointx, _pointy, _dirx, _diry)
      # if I am going in the right direction, then x would probably be the next point,
      # otherwise, I still move in the correct direction
       x = abs(x)
      # initial bracket
      ax= tol
      bx= 2*x
      __setbracket!(zb, ax, x, bx) # update the bracket
    end
    # solve
    findraymin(zb) # find the minimum


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

    #########################################################################################
    # update indices and normal direction
    ########################################################################################
    # initial normal assuming descending outward normal to a point on the Earth surface
    # need to compute the normal vector using cosd and sind because the θc might be clamped
    cosθc = cosd(θc)
    sinθc = sind(θc)
    N₀ = a/sqrt(1-e²*sinθc*sinθc) # normal to the Earth surface
    # explicityly used when clamped to find new t,h
    earthx = N₀*cosθc
    earthy = c²*N₀*sinθc
    ########################################################################################
    normalx= earthx/a²
    normaly= earthy/b²
    norm_hyp= hypot(normalx, normaly) # normalize the normal vector
    normalx /= norm_hyp
    normaly /= norm_hyp
    ########################################################################################
    if iter == 1
      if !isdescending
         return false # if it is the first iteration and not descending, return false
      end
      idx_j = 1
      idx_i = mod1(searchsortedlast(knots_θ, θc),N)
      # we know the ray is descending already
    elseif θc == θmin || θc == θmax
      # Find the correct value of t by using the intersection between two lines
      islevel = false

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

    elseif (abs(fmin)< ACCEPTABLE_TOLERANCE) ||
      hk< hmin ||
      hk > hmax
      if isdescending
        idx_j += 1 # set the next index
      else
        idx_j -= 1 # set the previous index
        normalx,normaly = -normalx, -normaly # reverse the normal direction
      end
    end
    ########################################################################################


    # update the point position
    rx +=  len_t * dirx
    ry +=  len_t * diry
    ########################################################################################
    # snells law (set to freespace if outside the grid)
    n_transmitted = (1<=idx_j<=M) ?  refractive_grid[idx_i,idx_j] : FREESPACE # get the refractive index at the current
    # update the direction

    pre_dirx = dirx # store the previous direction for logging
    pre_diry = diry # store the previous direction for logging
    dirx,diry = snellslaw(normalx, normaly, dirx, diry, n_i, n_transmitted)

    ########################################################################################
    # UPDATE FOR NEW ITERATION
    # update the refractive index for the next iteration
    n_i = n_transmitted # update the refractive index for the next iteration
    # h in descending order
    hmax = if 1<= idx_j <= M
      knots_h[idx_j]
    elseif idx_j < 1
      knots_h[1] # if idx_j is out of bounds, use the first value
    else
      knots_h[M] # if idx_j is out of bounds, use the last value
    end

    hmin = knots_h[idx_j+1] # get the hmax and hmin from the knots
    # note: knots_θ has an extra element to account for the last angle
    # theta in ascending order
    θmin = knots_θ[idx_i]
    θmax = knots_θ[mod1(idx_i+1,N)] # knots_θ has an extra element to account for the last angle
    # update for next iteration

    __update_ray!(zb,len_t, rx, ry, dirx, diry,
      hmin, hmax, θmin, θmax,
      idx_i, idx_j,n_i,
      isdescending, islevel) # update the ray
    # check if still in the atmosphere otherwise break the loop
    if ((idx_j < 1 || idx_j >= M) || len_t<tol)
      __update_ray!(zb,len_t, rx, ry, dirx, diry,
        hmin, hmax, θmin, θmax,
        idx_i, idx_j, n_i,
        isdescending, islevel) # update the ray

      return false
    end
    return true
end

#not inplace solve function
function solvenext(iter::Int,zb::Z,knots_h::V,knots_θ::V,
    refractive_grid::MA
) where {F,T<:AbstractFloat,V<:AbstractVector{T},
  MA<:AbstractMatrix{T},Z<:Zbrent{F,T}
}
    zb1=deepcopy(zb) # create a copy of the zbrent structure
    flag = solvenext!(iter,zb1,knots_h,knots_θ,refractive_grid) # solve the next step
    return zb1,flag # return the new zbrent structure
end



@inline function __solve!(res::A, zb::Z,
  pointx::T,pointy::T,
  directionx::T,directiony::T,
  knots_h::Vh,knots_θ::Vθ,
  refractive_grid::MA
) where {F,T<:AbstractFloat,R<:AbstractResult{T},
  A<:AbstractVector{R},Vh<:AbstractVector{T},Vθ<:AbstractVector{T},
  MA<:AbstractMatrix{T},Z<:Zbrent{F,T}
}

  N= length(knots_θ)-1
  M= length(knots_h)-1
  niters= length(res)-1 # number of iterations to perform
  # Initialize the ray structure
  hmin = knots_h[1]
  n_i = FREESPACE

  solveinit!(zb,n_i, pointx, pointy, directionx, directiony,hmin)
  __update_results!(res, 0, zb) # update the results
  @inbounds for iter in 1:niters
    # find the minimium in the bracket

    ok = solvenext!(iter,zb,view(knots_h,:),view(knots_θ,:),view(refractive_grid,:,:)) # solve the next step

        # update the results
    ########################################################################################
    __update_results!(res, iter, zb) # update the results
    ########################################################################################

    if !ok
      break # if the ray is not in the atmosphere, break the loop
    end
  end
end
