
# create bowring approximation for geodetic coordinates
begin
  ##########################################################################################################################################
  local name=:bowring
  local return_angle = Expr[:(mod1(atand(Num,Den*ang_sign),360))]
  local return_altitude = Expr[:(tau2=tau*tau), :((W+Z*tau-c5*sqrt(1+c2*tau*tau))/sqrt(1+tau*tau))]
  local return_altitudeangle = Expr[ :(tau2=tau*tau), :((W+Z*tau-c5*sqrt(1+c2*tau*tau))/sqrt(1+tau*tau),mod1(atand(Num,Den*ang_sign),360))]
  ##########################################################################################################################################
  local RETURNTHIS = (return_angle, return_altitude, return_altitudeangle)
  for (what,out,ret,safe) in zip(RETURNWHAT,OUTTYPE,RETURNTHIS,SAFERETURN)
    fun = Symbol("ray2_",name,what)
    @inline @eval  function $fun(W::T,Z::T)::$(out) where T<:AbstractFloat
        $(Expr(:block, safe...))
        ang_sign= sign(W)
        c2 =COMPLEMENTECCENTRICITY²(T)
        c1 = 1/c2
        c3 = c1*c1
        c4 = MAJORAXIS(T) * (1 - c2)
        c5 = MAJORAXIS(T)

        W= abs(W)
        W2   = W*W
        Z2= Z*Z

        K    = W2+c1*Z2
        L    = c4/(K*sqrt(K))
        Num  = Z + c3*Z2*Z*L
        Den  = W - W2*W*L
        tau  = Num / Den
        $(Expr(:block, ret...))
    end
    @inline @eval $fun(args) = $fun(args...)
  end
