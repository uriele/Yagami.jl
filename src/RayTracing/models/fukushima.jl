
# create fukushima approximation for geodetic coordinates
  name=:fukushima
  return_angle = Expr[:(mod1(atand(s1,cc*ang_sign),360))]
  return_altitude = Expr[:(s12=s1*s1),:(cc2=cc*cc),
  :((W  *cc+Z*s1-   c6*sqrt(c2*s12   +cc2   ))/sqrt(s12+cc2))]
  return_altitudeangle = Expr[:(s12=s1*s1),:(cc2=cc*cc),
  :(((W  *cc+Z*s1-   c6*sqrt(c2*s12   +cc2   ))/sqrt(s12+cc2),mod1(atand(s1,cc*ang_sign),360)))]
  for (what,out,ret,safe) in zip((:_angle,:_altitude, :_altitudeangle),
    outtype,(return_angle,return_altitude,return_altitudeangle),safe_return)
    fun = Symbol("ray2_",name,what)
    @inline @eval  function $fun(W::T,Z::T)::$(out) where T<:AbstractFloat
        ## safeguard
        $(Expr(:block, safe...))

        ## Fukushima approximation for geodetic coordinates
        ang_sign= sign(W)
        c1=1/MAJORAXIS(T)
        c2=COMPLEMENTECCENTRICITY²(T)
        c3=1-c2
        c4=sqrt(c2)
        c5=1.5*c3*c3
        c6=MAJORAXIS()
        W = abs(W)
        Z = Z
        s0=c1*Z
        Wn=c1*W
        c0=c4*Wn
        c02=c0*c0
        s02=s0*s0
        a02=c02+s02
        a0=sqrt(a02)
        a03=a02*a0
        f0=Wn*a03   -c3*c02   *c0
        b0=c5*s02   *c02   *Wn*(a0-c4)
        s1=(c4*s0*a03+c3*s02   *s0)*f0-b0*s0
        cc=c4*(f0*f0-b0*c0)
        $(Expr(:block, ret...))
    end
    @inline @eval $fun(args) = $fun(args...)
end
