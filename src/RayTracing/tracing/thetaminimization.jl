using ..YagamiCore:difference_angles
#initialize the solver
@inline function solveinitθ!(zb::ZB,n_i::T,pointx::T,pointy::T,
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
  # This must be set to avoid that the initial angle would be on the opposite side of the
  # Earth surface, forcing the minimum to the be on the wrong side of the Earth.
  θ₀= atand(pointy,pointx) #initial angle approximated to the position on a circle.
  __setbracket!(zb, θ₀,0.0 , 0.0) # reset the bracket
end

@inline function newton_step!(zb::ZB) where {F,T<:AbstractFloat,ZB<:Zbrent{F,T}}
    ###########################################################
    # earth model parameters
    ###########################################################
    a_orig = MAJORAXIS(T)
    c² = COMPLECCENTRICITY²(T)
    e² = ECCENTRICITY²(T)

    tol = zb.tol
    rx = getpointx(zb)/a_orig
    ry = getpointy(zb)/a_orig
    dirx = getdirectionx(zb)
    diry = getdirectiony(zb)

    isdescending = getdescending(zb)
    hmin = gethmin(zb)/a_orig
    htest = hmin

    θ = zb.a # initial angle, slightly offset to avoid zero division

    a= one(T) # major axis
    dist²= zero(T) # distance squared
    t= zero(T) # distance along the ray
    for i in 0:zb.itermax
        cosθ = cosd(θ)
        sinθ = sind(θ)
        sin2θ_half = sinθ*cosθ
        cosθ² = cosθ*cosθ
        sinθ² = sinθ*sinθ
        cos2θ = cosθ² - sinθ²

        R₀ = sqrt(1-e²*sinθ²)
        invR₀³ = 1/(R₀^3)
        invR₀⁵ = 1/(R₀^5)
        N₀ = a/R₀
        dN₀ = e²*a*sin2θ_half*invR₀³
        d²N₀ = e²*a*(cos2θ*invR₀³ - 3*e²*sin2θ_half^2*invR₀⁵)

        Fx0 = N₀*cosθ
        Fy0 = c²*N₀*sinθ

        Fx1 = cosθ
        Fy1 = sinθ

        Fx = Fx0 + htest*Fx1
        Fy = Fy0 + htest*Fy1

        dFxdθ = dN₀*cosθ - sinθ*(N₀+htest)
        dFydθ = c²*dN₀*sinθ + cosθ*(c²*N₀+htest)

        t = (Fx-rx)*dirx + (Fy-ry)*diry

        Px = rx + t*dirx
        Py = ry + t*diry

        FPx = Fx - Px
        FPy = Fy - Py

        dist² = FPx^2 + FPy^2

        dt_dθ = (dFxdθ*dirx + dFydθ*diry)
        d²t_dθ² = (
                  (d²N₀*cosθ - 2*dN₀*sinθ - (N₀+htest)*cosθ)*dirx +
                  (c²*d²N₀*sinθ + 2*c²*dN₀*cosθ - (c²*N₀+htest)*sinθ)*diry
                  )

        d_distx = dFxdθ - dt_dθ*dirx
        d_disty = dFydθ - dt_dθ*diry

        gstep = d_distx*FPx + d_disty*FPy

        d²Fxdθ² = d²N₀*cosθ - 2*dN₀*sinθ - (N₀+htest)*cosθ
        d²Fydθ² = c²*d²N₀*sinθ + 2*c²*dN₀*cosθ - (c²*N₀+htest)*sinθ

        hstep = d_distx^2 + d_disty^2 + FPx*(d²Fxdθ² - dirx*d²t_dθ²) +
                 FPy*(d²Fydθ² - diry*d²t_dθ²)

        if abs(hstep) < eps(T) || isnan(hstep) || isnan(gstep)
            @warn "Hessian ill-conditioned or NaN appeared at iteration $i with d²fdθ²: $d²fdθ², gstep: $gstep"
            break # Exit if Hessian is ill-conditioned or NaN appears
        end

        pstep = gstep / hstep

        θ = mod(θ - pstep, 360)

        if sqrt(dist²) < tol
            zb.x = t*a_orig
            zb.a = mod(θ,360) # update the angle in the Zbrent structure
            zb.fx = sqrt(dist²)*a_orig
            return
        end
    end

    zb.x = t*a_orig
    zb.a = mod(θ,360)
    zb.fx = sqrt(dist²)*a_orig
end


@inline function newton_zero!(zb::ZB) where {F,T<:AbstractFloat,ZB<:Zbrent{F,T}}
    ###########################################################
    # earth model parameters
    ###########################################################
    a = MAJORAXIS(T)
    c² = COMPLECCENTRICITY²(T)
    e² = ECCENTRICITY²(T)

    tol = zb.tol
    rx = getpointx(zb)
    ry = getpointy(zb)
    dirx = getdirectionx(zb)
    diry = getdirectiony(zb)

    t = 1e-5 # small offset to avoid zero division

    hmax = gethmax(zb)
    htest = hmax

    θ = zb.a # initial angle, slightly offset to avoid zero division

    _,θ1=__innerfunc(zb,rx+t*dirx, ry+t*diry) # compute the inner function at the initial point

    θ+=sign(difference_angles(θ1,θ)) # adjust the initial angle to the inner function value

    dist = zero(T) # distance squared
    t= zero(T) # distance along the ray

    for i in 0:zb.itermax
      Fx1=cosθ = cosd(θ)
      Fy1=sinθ = sind(θ)
      sin2θ_half = sinθ*cosθ
      cosθ² = cosθ*cosθ
      sinθ² = sinθ*sinθ

      R₀ = sqrt(1-e²*sinθ²)
      invR₀³ = 1/(R₀^3)

      N₀ = a/R₀
      dN₀ = e²*a*sin2θ_half*invR₀³

      Fx0 = N₀*cosθ
      Fy0 = c²*N₀*sinθ

      Δx = Fx0 - rx
      Δy = Fy0 - ry

      num=(Δy*dirx - Δx*diry)

      det = -dirx*Fy1 + diry*Fx1

      det= abs(det)<eps(T) ? eps(T)*sign(det) : det

      h = num/det
      t = (Δy*Fx1 - Δx*Fy1)/det
      dist = htest-h

      if abs(dist) < tol
        zb.x = t
        zb.a = mod(θ,360) # update the angle in the Zbrent structure
        zb.fx = htest - h
        return
      end

      dFx0dθ = dN₀*cosθ-sinθ*N₀
      dFy0dθ = c²*(dN₀*sinθ +N₀* cosθ)
      dFx1dθ = -sinθ
      dFy1dθ = cosθ
      dnum= dFy0dθ*dirx-dFx0dθ*diry
      dden= (dirx*dFy1dθ - diry*dFx1dθ)/det
      ddistdθ = -1/det*(dnum+dden)
      dθ= dist/ddistdθ

      θ = mod(θ - dθ, 360)


    end
    zb.x = t
    zb.a = mod(θ,360) # update the angle in the Zbrent structure
    zb.fx = dist


end

@inline function __update_rayθ!(zb::Z,len_t::T,θc::T,hk::T,pointx::T,pointy::T,
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
  __setbracket!(zb, θc, len_t, hk) # reset the bracket
end

#inplace solve function
function solvenextθ!(iter::Int,zb::Z,
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

    # clamp the angle to the range [θmin, θmax]
    if isdescending
      newton_step!(zb) # update the angle in the Zbrent structure
    else
      newton_zero!(zb) # update the angle in the Zbrent structure
    end

    len_t = zb.x # get the length of the ray
    fmin  = zb.fx
    # trial
    _x = rx+len_t * dirx
    _y = ry+len_t * diry
    hk,θk = __getinnerf(zb, _x, _y) # get the height and angle at the point
    # check if it is in the range (θmin, θmax)

    θc = clampangle(θk, θmin, θmax)

    zb.a = θc # update the angle in the Zbrent structure
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
      hk<= hmin+ACCEPTABLE_TOLERANCE ||
      hk >= hmax+ACCEPTABLE_TOLERANCE
      if isdescending
        idx_j += 1 # set the next index
      else
        idx_j -= 1 # set the previous index
        normalx,normaly = -normalx, -normaly # reverse the normal direction
      end
    else
      isdescending = false # reverse the direction
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

    __update_rayθ!(zb,len_t,θc,hk, rx, ry, dirx, diry,
      hmin, hmax, θmin, θmax,
      idx_i, idx_j,n_i,
      isdescending, islevel) # update the ray
    # check if still in the atmosphere otherwise break the loop
    if ((idx_j < 1 || idx_j >= M) || len_t<tol)
      __update_rayθ!(zb,len_t,θc,hk, rx, ry, dirx, diry,
        hmin, hmax, θmin, θmax,
        idx_i, idx_j, n_i,
        isdescending, islevel) # update the ray

      return false
    end
    return true
end




@inline function __solveθ!(res::A, zb::Z,
  pointx::T,pointy::T,
  directionx::T,directiony::T,
  knots_h::Vh,knots_θ::Vθ,
  refractive_grid::MA
) where {F,T<:AbstractFloat,R<:AbstractResult{T},
  A<:AbstractVector{R},Vh<:AbstractVector{T},Vθ<:AbstractVector{T},
  MA<:AbstractMatrix{T},Z<:Zbrent{F,T}
}

  M= length(knots_h)-1
  niters= length(res)-1 # number of iterations to perform
  # Initialize the ray structure
  hmin = knots_h[1]
  n_i = FREESPACE

  solveinitθ!(zb,n_i, pointx, pointy, directionx, directiony,hmin)
  __update_results!(res, 0, zb) # update the results
  @inbounds for iter in 1:niters
    # find the minimium in the bracket

    ok = solvenextθ!(iter,zb,view(knots_h,:),view(knots_θ,:),view(refractive_grid,:,:)) # solve the next step

        # update the results
    ########################################################################################
    __update_results!(res, iter, zb) # update the results
    ########################################################################################

    if !ok
      break # if the ray is not in the atmosphere, break the loop
    end
  end
end




"""
    $SIGNATURES
Perform ray tracing calculations for a given problem and store the results in the provided results matrix. This function iterates over the points and directions specified in the problem, solving for the ray path using the Zbrent method.
It updates the results matrix with the ray tracing results for each point and direction. It also writes a log of the ray tracing results to the logger associated with the problem.

# Arguments:
- `results::RR`: The results matrix where the ray tracing results will be stored.
- `problem::RP`: The ray tracing problem containing the points, directions, refractive index grid, and other necessary data.
- `itermax::Int`: The maximum number of iterations to perform for each ray tracing calculation (default is 100).
- `tol::T`: The tolerance for the ray tracing calculations (default is 1e-10).
"""
function raytracingθ!(
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
  knots_h = view(problem.atmosphere.temperature.knots_h,:)
  knots_θ = view(problem.atmosphere.temperature.knots_θ,1:Ni)
  refractive_grid = problem.refractive

  model = problem.earthmodel

  logger = getlogger(problem)
  nlos = problem.nlos
  nscans = problem.nscans

  zb = Zbrent(T,DistanceFunc(earthmodelfunction(model)),itermax,tol)
  @inbounds for i in eachindex(pointsx)
    @inline __solveθ!(view(results,:,i),zb, pointsx[i], pointsy[i], directionsx[i], directionsy[i],knots_h, knots_θ, refractive_grid)
  end

  write_tracing_log(logger, results,refractive_grid)
end


"""
    $SIGNATURES
Perform parallel ray tracing calculations for a given problem and store the results in the provided results matrix. This function uses multithreading to perform ray tracing calculations for each point and direction specified in the problem.
It updates the results matrix with the ray tracing results for each point and direction. It also writes a log of the ray tracing results to the logger associated with the problem.
# Arguments:
- `results::RR`: The results matrix where the ray tracing results will be stored.
- `problem::RP`: The ray tracing problem containing the points, directions, refractive index grid, and other necessary data.
- `itermax::Int`: The maximum number of iterations to perform for each ray tracing calculation (default is 100).
- `tol::T`: The tolerance for the ray tracing calculations (default is 1e-10).
"""
function raytracingθ_parallel!(
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
  knots_h = view(problem.atmosphere.temperature.knots_h,:)
  knots_θ = view(problem.atmosphere.temperature.knots_θ,1:Ni)
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
    @inbounds @inline __solveθ!(view(results,:,i),zb_buffer[Threads.threadid()], pointsx[i], pointsy[i], directionsx[i], directionsy[i],knots_h, knots_θ, refractive_grid)
  end

  write_tracing_log(logger, results,refractive_grid)
end
