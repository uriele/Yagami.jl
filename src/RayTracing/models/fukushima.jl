
# create fukushima approximation for geodetic coordinates
begin
  local name=:fukushima

   ##########################################################################################################################################
  local return_angle = Expr[ FUKUSHIMA_RETURN_ANGLE]
  local return_altitude = Expr[ FUKUSHIMA_RETURN_ALTITUDE...]
  local return_altitudeangle = Expr[FUKUSHIMA_RETURN_ALTITUDEANGLE...]
  ##########################################################################################################################################
  local RETURNTHIS = (return_angle, return_altitude, return_altitudeangle)

  for (what,out,ret,safe) in zip(RETURNWHAT,OUTTYPE,RETURNTHIS,SAFERETURN)


    # inline function for sensitive calculations
    inlinefun = Symbol("__",what,"_",name)

    @inline @eval function $inlinefun(W::T,Z::T,_MAJORAXIS::T,_MINORAXIS::T,_COMPLECCENTRICITY²::T) where T<:AbstractFloat
        $(Expr(:block, safe...))
        $(Expr(:block, FUKUSHIMA_MAIN_BODY...))
        $(Expr(:block, ret...))
    end


    # main function for the fukushima model
    fun = Symbol("ray2",what,"_",name)
    @inline @eval  function $fun(W::T,Z::T)::$(out) where T<:AbstractFloat
        _MAJORAXIS = MAJORAXIS(T)
        _MINORAXIS = MINORAXIS(T)
        _COMPLECCENTRICITY² = COMPLECCENTRICITY²(T)
        $inlinefun(W,Z,_MAJORAXIS,_MINORAXIS,_COMPLECCENTRICITY²)
    end


    # verbose version for plotting and debugging
    fun_verbose = Symbol("ray2",what,"_",name,"_verbose")
    @eval  function $fun_verbose(W::T,Z::T,a::T=WGS84MAJORAXIS,b::T=WGS84MINORAXIS)::$(out) where T<:AbstractFloat
        _COMPLECCENTRICITY²=b^2/a^2
        $inlinefun(W,Z,a,b,_COMPLECCENTRICITY²)
    end

    @eval $fun(args) = $fun(args...)
    @eval $fun_verbose(args) = $fun_verbose(args...)
  end

end
