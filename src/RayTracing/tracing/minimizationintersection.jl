const ACCEPTABLE_TOLERANCE::Float64 =1e-5
const FREESPACE::Float64 = 1.0 # free space index
struct NoFunc end
# pass traits of the Earth Approximation
EarthApproximation(x) = EarthApproximation(typeof(x)) # convert to type
EarthApproximation(::Type)= NoFunc() # no type conversion
EarthApproximation(::T) where {T<:Type{<:Fukushima}} = Fukushima
EarthApproximation(::T) where {T<:Type{<:Bowring}} = Bowring

@inline earthmodelfunction(x)=earthmodelfunction(EarthApproximation(x))
@inline earthmodelfunction(::NoFunc) = error("No Earth Approximation specified. Use Fukushima or Bowring.")
@inline earthmodelfunction(::Type{<:Fukushima}) = ___altitudeangle_fukushima

@inline earthmodelfunction(::Type{<:Bowring}) = ___altitudeangle_bowring


@inline __innerfunc(z::Zbrent, x::T,y::T) where T = getfield(z.f, :__f)(x,y,MAJORAXIS(T),MINORAXIS(T),COMPLECCENTRICITY²(T))


function raytracing!(results::RR,
  pointsx::V,pointsy::V,directionsx::V,directionsy::V,
  knots_h::Vh,knots_θ::Vθ,refractive_grid::M,model::EA=Fukushima,itermax::Int=100,tol=1e-10) where {T<:AbstractFloat,R<:ResultRay{T},
  EA<:EarthApproximation,
  RR<:AbstractVector{R},V<:AbstractVector{T}, Vh<:AbstractVector{T},
  Vθ<:AbstractVector{T}, M<:AbstractMatrix{T}}

  zb = Zbrent(T,DistanceFunc(earthmodelfunction(model)),itermax,tol)

  @inbounds for i in eachindex(pointsx)
    @inline __solve!(res[i],zb, pointsx[i], pointsy[i], directionsx[i], directionsy[i],view(knots_h,:), view(knots_θ,:), view(refractive_grid,:,:))
  end
end





@inline function __solve!(res::ResultRay, zb::Zbrent{F,T},
  pointx::T,pointy::T,directionx::T,directiony::T,
  knots_h::AbstractVector{T},knots_θ::AbstractVector{T},
  refractive_grid::AbstractMatrix{T}) where {F,T<:AbstractFloat}

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

  __sethlims!(zb, hmin, hmax)
  __setθlims!(zb, θmin, θmax)
  ###########################
  __setdescending!(zb) # reset to descending
  for i in 1:niters
    bracket(zb)
    if zb.b< zb.tol
      __setdescending!(zb,false) # reset to descending
      bracket(zb) # recalculate the bracket
    end
    findmin(zb) # find the minimum
    hk,θk = __getinnerf(zb.f, zb.pointx, zb.pointy)
    # check if it is in the range (θmin, θmax)
    θc = clampangle(θk, θmin, θmax)

    # normal direction with respect to the intersection point
    # outward normal for the descending ray
    normalx = cosd(θc)
    normaly = sind(θc)

    if iter == 1
      getdescenting(zb) || break # if not descending, break the loop it means that it never entered the atmosphere
      __setj!(zb, 1) # set the index for the first iteration
      # find the radial element
      # find the last element that is less than or equal to θc
      i = searchsortedlast(view(knots_θ,1:N), θc)
      __seti!(zb, i) # set the index
      # we know the ray is descending already
    end
    if θc == θmin || θc == θmax
      # Find the correct value of t by using the intersection between two lines
      No = MAJORAXIS(T)/sqrt(1-ECCENTRICITY²(T)*normaly*normaly) # normal to the Earth surface
      earthx = N₀*normalx
      earthy = COMPLECCENTRICITY²(T)*N₀*normaly

      t,h = intersectionrayray(pointx, pointy, directionx, directiony, earthx, earthy, normalx, normaly)

      # bracket contains z.a,z.x,z.b,z.f(x)
      __setbracket!(zb, t, t, t) # set the bracket to the intersection point

      # important: this is the direction for θc=θmax
      normalx,normaly = normaly, -normalx # swap the normal direction
      __setislevel(zb, false) # set to islevel
      isθmax = θc == θmax
      if isθmax
        __seti!(zb,mod(i+1,N)) # set the next index
      else
        __seti!(zb,mod(i-1,N)) # set the previous index
        normalx,normaly = -normalx, -normaly # reverse the normal direction
      end

    elseif zb.fx< ACCEPTABLE_TOLERANCE
      __setislevel(zb, true) # set to islevel
      if getdescending(zb)
        __setj!(zb,i+1) # set the next index
      else
        __setj!(zb,i-1) # set the previous index
        normalx,normaly = -normalx, -normaly # reverse the normal direction
      end
    end
    # if the ray stays inside the wedge, we have a minimum
    # so the normal is not important since n_i==n_t and the
    # Snell's law is not applied

    # check if still in the atmosphere otherwise break the loop
    if geti(zb) < 1 || geti(zb) > M
      break # if the index is out of bounds, break the loop
    end

    # update the point position
    pointx += zb.x * zb.directionx
    pointy += zb.x * zb.directiony

    n_transmitted = refractive_grid[geti(zb), getj(zb)] # get the refractive index at the current

    # update the direction
    directionx,directiony = snellslaw(normalx, normaly, directionx, directiony, n_i, n_transmitted)

    hmax,hmin = knots_h[getj(zb)],knots_h[geti(zb)+1] # get the hmax and hmin from the knots
    # note: knots_θ has an extra element to account for the last angle
    θmax,θmin = knots_θ[getj(zb)],knots_θ[geti(zb)+1] # get the θmax and θmin from the knots
    # update for next iteration
    __setpoint!(zb, pointx, pointy) # set the new point position
    __setdirection!(zb, directionx, directiony) # set the new direction
    __sethlims!(zb, hmin, hmax) # set the new h limits
    __setθlims!(zb, θmin, θmax) # set the new θ limits
  end

end
