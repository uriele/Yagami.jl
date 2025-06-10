const OUTTYPE = [:(T),:(T),:(Tuple{T,T})]
const SAFEANGLE = Expr[ :(W==0 && return Z==0 ? 360 : (Z>0 ? 90 : 270)),
  :(Z == 0 && return W>0 ? 360 : 180)]
const SAFEALTITUDE = Expr[ :(W==0 && return Z==0 ? 0 : abs(Z)-MINORAXIS(T)),
  :(Z==0 && return W==0 ? 0 : abs(W)-MAJORAXIS(T))]
const SAFEALTITUDEANGLE = Expr[
    :(W==0 && return Z==0 ? (0,0) : (abs(Z)-MINORAXIS(T),Z>0 ? 90 : 270)),
    :(Z==0 && return W==0 ? (0,0) : (abs(W)-MAJORAXIS(T),W>0 ? 360 : 180))]
const SAFERETURN = (SAFEANGLE,SAFEALTITUDE,SAFEALTITUDEANGLE)

const RETURNWHAT = (:_angle,:_altitude, :_altitudeangle)
############################################################################################
# Constants for Bowring approximation
############################################################################################safe_returns_angle = Expr[ :(W==0 && return Z==0 ? 360 : (Z>0 ? 90 : 270)),
