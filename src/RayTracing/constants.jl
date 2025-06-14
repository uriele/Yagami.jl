const SMALLSHIFT= 1e-6
############################################################################################
const OUTTYPE = [:(T),:(T),:(Tuple{T,T})]
const SAFEANGLE = Expr[ :(W==0 && return Z==0 ? 360 : (Z>0 ? 90 : 270)),
  :(Z == 0 && return W>0 ? 360 : 180)]
const SAFEALTITUDE = Expr[ :(W==0 && return Z==0 ? 0 : abs(Z)-_MINORAXIS),
  :(Z==0 && return W==0 ? 0 : abs(W)-_MAJORAXIS)]
const SAFEALTITUDEANGLE = Expr[
    :(W==0 && return Z==0 ? (0,0) : (abs(Z)-_MINORAXIS,Z>0 ? 90 : 270)),
    :(Z==0 && return W==0 ? (0,0) : (abs(W)-_MAJORAXIS,W>0 ? 360 : 180))]
const SAFERETURN = (SAFEANGLE,SAFEALTITUDE,SAFEALTITUDEANGLE)

const RETURNWHAT = (:_angle,:_altitude, :_altitudeangle)
############################################################################################
# Constants for Bowring approximation
############################################################################################safe_returns_angle = Expr[ :(W==0 && return Z==0 ? 360 : (Z>0 ? 90 : 270)),
const BOWRING_MAIN_BODY = Expr[
        :(ang_sign= sign(W)),
        :(c2 =_COMPLECCENTRICITY²),
        :(c1 = 1/c2),
        :(c3 = c1*c1),
        :(c4 = _MAJORAXIS * (1 - c2)),
        :(c5 = _MAJORAXIS),

        :(W= abs(W)),
        :(W2   = W*W),
        :(Z2= Z*Z),

        :(K    = W2+c1*Z2),
        :(L    = c4/(K*sqrt(K))),
        :(Num  = Z + c3*Z2*Z*L),
        :(Den  = W - W2*W*L),
        :(tau  = Num / Den)
]

const BOWRING_RETURN_ANGLE    = :(mod(atand(Num,Den*ang_sign),360))

const BOWRING_RETURN_ALTITUDE = :((W+Z*tau-c5*sqrt(1+c2*tau*tau))/sqrt(1+tau*tau))

const BOWRING_RETURN_ALTITUDEANGLE = :(return ((W+Z*tau-c5*sqrt(1+c2*tau*tau))/sqrt(1+tau*tau)),mod1(atand(Num,Den*ang_sign),360))
############################################################################################
# Constants for Fukushima approximation
############################################################################################

const FUKUSHIMA_MAIN_BODY = Expr[
        :(ang_sign= sign(W)),
        :(c1=1/_MAJORAXIS),
        :(c2=_COMPLECCENTRICITY²),
        :(c3=1-c2),
        :(c4=sqrt(c2)),
        :(c5=1.5*c3*c3),
        :(c6=_MAJORAXIS),
        :(W = abs(W)),
        :(Z = Z),
        :(s0=c1*Z),
        :(Wn=c1*W),
        :(c0=c4*Wn),
        :(c02=c0*c0),
        :(s02=s0*s0),
        :(a02=c02+s02),
        :(a0=sqrt(a02)),
        :(a03=a02*a0),
        :(f0=Wn*a03   -c3*c02   *c0),
        :(b0=c5*s02   *c02   *Wn*(a0-c4)),
        :(s1=(c4*s0*a03+c3*s02   *s0)*f0-b0*s0),
        :(cc=c4*(f0*f0-b0*c0))
]

const FUKUSHIMA_RETURN_ANGLE     = :(mod(atand(s1,cc),360))
const FUKUSHIMA_RETURN_ALTITUDE  = :(s12=s1*s1),:(cc2=cc*cc),:((W  *cc+Z*s1-   c6*sqrt(c2*s12   +cc2   ))/sqrt(s12+cc2))
const FUKUSHIMA_RETURN_ALTITUDEANGLE = Expr[
    :(s12=s1*s1),
    :(cc2=cc*cc),
    :(return ((W  *cc+Z*s1-   c6*sqrt(c2*s12   +cc2   ))/sqrt(s12+cc2),mod1(atand(s1,cc*ang_sign),360)))]
############################################################################################

const EXISTINGMODELS = (:fukushima, :bowring)


# Constants for the Earth approximations
const EXPFUNC = (:getpoint,:getdirection, :getwedgeindex, :getlength, :getaltitude, :getazimuth)
