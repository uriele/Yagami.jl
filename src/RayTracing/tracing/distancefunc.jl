
mutable struct DistanceFunc{F,T<:AbstractFloat}
    # parameters needed for the zberent method
    __f::F
    pointx::T
    pointy::T
    directionx::T
    directiony::T
    hmin::T
    hmax::T
    # used to check if we are in bounding box,
    # happens afret minimization
    θmin::T
    θmax::T
    # store current i,j to pass easily the information about the wedge
    i ::Int
    j:: Int
    n:: T
    # use to compute the correct distance for the iteration and for the normal
    descending::Bool
    # find automatically if it hits a level or a radii to compute the normal
    islevel::Bool
    #
    @inline function DistanceFunc(f::F=___altitudeangle_fukushima, pointx::T=DFPOINTX, pointy::T=DFPOINTY,
                                     directionx::T=DFDIRECTIONX,directiony::T=DFDIRECTIONY,
                                     hmin::T=DFHMIN, hmax::T=DFHMAX,
                                     θmin::T=DFTHMIN, θmax::T=DFTHMAX,n::T=FREESPACE) where {F,T<:AbstractFloat}

      # Normalize direction vector
      direction_norm = hypot(directionx, directiony)
      directionx /= direction_norm
      directiony /= direction_norm

      return new{F,T}(f, pointx, pointy, directionx, directiony, hmin, hmax, θmin, θmax,
                      DFI, DFJ,n, DFDESCENDING, DFISLVL)
    end
end

Base.show(io::IO, d::DistanceFunc{F,T}) where {F,T} = print(io, "DistanceFunc{T}")
Base.propertynames(::DistanceFunc) = DISTFUNCPROP

@inline function Base.getproperty(df::DistanceFunc, sym::Symbol)
  if sym ∈ DISTFUNCPROP
      return getfield(df, sym)
  end
end

# cannot set properties of DistanceFunc directly to avoid issues
@inline function Base.setproperty!(::DistanceFunc, ::Symbol, value)
  return nothing
end
############################################################################################
# Getters and Setters for DistanceFunc
for (propr,defval,proptype) in zip(DISTFUNCPROP, DFDISTFUNCPROP, DISTTYPEPROP)
    getfunc = Symbol("get", propr)

    docstr_get = """
    `$getfunc(df::DistanceFunc)`
    Return the `:$propr` field of the `DistanceFunc`.
    """
    @eval begin
        @doc $docstr_get
        @inline function $getfunc(df::DistanceFunc)
            getfield(df, $(QuoteNode(propr)))
        end

        @forward Zbrent.f $getfunc
    end

    setfunc = Symbol("__set", propr, "!")
    docstr_set = """
    `$setfunc(df::DistanceFunc, value::$proptype=$defval)`
    Set the `:$propr` field of the `DistanceFunc`.
    """
    @eval begin
        @doc $docstr_set
        @inline function $setfunc(df::DistanceFunc, value::$proptype = $defval) where {$proptype}
            setfield!(df, $(QuoteNode(propr)), value)
        end
        @forward Zbrent.f $setfunc
    end
end
# easier larger functions
getpoint(df)    = getpointx(df), getpointy(df)
getdirection(df)= getdirectionx(df), getdirectiony(df)
gethlims(df) = gethmin(df), gethmax(df)
getθlims(df) = getθmin(df), getθmax(df)

# easier setter

@inline function __setpoint!(df::DistanceFunc, pointx::T, pointy::T) where {T<:AbstractFloat}
  __setpointx!(df, pointx)
  __setpointy!(df,  pointy)
end

@inline function __setdirection!(df::DistanceFunc, directionx::T, directiony::T) where {T<:AbstractFloat}
__setdirectionx!(df, directionx)
__setdirectiony!(df, directiony)

end
@inline function __sethlims!(df::DistanceFunc, hmin::T, hmax::T) where {T<:AbstractFloat}
__sethmin!(df, hmin)
__sethmax!(df, hmax)
end

@inline function __setθlims!(df::DistanceFunc, θmin::T, θmax::T) where {T<:AbstractFloat}
__setθmin!(df, θmin)
__setθmax!(df, θmax)
end

@inline __getinnerf(df::DistanceFunc{F,T}, x::T, y::T) where {F,T<:AbstractFloat} =  getfield(df, :__f)(x, y, MAJORAXIS(T), MINORAXIS(T), COMPLECCENTRICITY²(T))
# Getters and Setters for DistanceFunc properties
@forward Zbrent.f getpoint, getdirection, gethlims, getθlims
@forward Zbrent.f __setpoint!, __setdirection!, __sethlims!, __setθlims!
@forward Zbrent.f __getinnerf


function (df::DistanceFunc{F,T})(t::T) where {F,T<:AbstractFloat}
  isDescending = getdescending(df)

  pointx = getpointx(df)
  pointy = getpointy(df)
  x   =   pointx + t * getdirectionx(df)
  y   =   pointy + t * getdirectiony(df)
  h,_ = __getinnerf(df,x ,y)

  htest= if isDescending
    gethmin(df)
  else
    gethmax(df)
  end
  dh=abs(htest-h)

end
