############################################################################################
# Constants for the DistanceFunctions
############################################################################################
# Default values for the distance function parameters
const DFI::Int = 0
const DFJ::Int = 0
const DFISLVL::Bool = true
const DFDESCENDING::Bool = true
const DFPOINTX::Float64 = 0.00001
const DFPOINTY::Float64 = 0.00001
const DFDIRECTIONX::Float64 = 0.0
const DFDIRECTIONY::Float64 = 1.0
const DFHMIN::Float64 = 0.0
const DFHMAX::Float64 = 120.0
const DFTHMIN::Float64 = -Inf
const DFTHMAX::Float64 = Inf
const DFFREESPACE::Float64 = 1.0
# note: f not directly accessible
const DISTFUNCPROP  = (:pointx,   :pointy,  :directionx, :directiony,  :hmin,  :hmax,   :θmin,   :θmax, :i, :j, :n, :descending, :islevel)
const DFDISTFUNCPROP      = (DFPOINTX, DFPOINTY, DFDIRECTIONX, DFDIRECTIONY,DFHMIN, DFHMAX, DFTHMIN, DFTHMAX,DFI,DFJ,DFFREESPACE,DFDESCENDING,  DFISLVL)
const DISTTYPEPROP        = (:T, :T, :T, :T, :T, :T, :T, :T, :Int, :Int,:T, :Bool, :Bool)

############################################################################################
const SMALLSHIFT= 1e-6
############################################################################################
const OUTTYPE = [:(T),:(T),:(Tuple{T,T})]
const SAFEANGLE = Expr[ :(W==0 && return Z==0 ? T(360) : (Z>0 ? T(90) : T(270))),
  :(Z == 0 && return W>0 ? 360 : 180)]
const SAFEALTITUDE = Expr[ :(W==0 && return Z==0 ? zero(T) : abs(Z)-_MINORAXIS),
  :(Z==0 && return W==0 ? zero(T) : abs(W)-_MAJORAXIS)]
const SAFEALTITUDEANGLE = Expr[
    :(W==0 && return Z==0 ? (zero(T),zero(T)) : (abs(Z)-_MINORAXIS,Z>0 ? T(90) : T(270))),
    :(Z==0 && return W==0 ? (zero(T),zero(T)) : (abs(W)-_MAJORAXIS,W>0 ? T(360) : T(180)))]
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
############################################################################################
# Constants for the PROBLEMRAY data structure
############################################################################################

const PROBLPROPDIRECT =(:filename, :meantype, :model, :earthmodel,
  :refractive, :pointsx, :pointsy, :atmosphere,:directionsx, :directionsy, :tangent_h, :tangent_θ,
  :nscans, :nlos)

const PROBLATMPROP = (:temperature,:pressure, :humidity, :co2ppm, :wavelength)
const PROBLKNOTS = (:knots_h, :knots_θ)
const PROBLPROPTOTAL= vcat(PROBLPROPDIRECT..., PROBLATMPROP..., PROBLKNOTS...)
############################################################################################
# Constants for the RayTracingProblem
#############################################################################################

const ACCEPTABLE_TOLERANCE::Float64 =1e-5
const FREESPACE::Float64 = 1.0 # free space index

############################################################################################
# Constants to read the Geofit format files
############################################################################################

# names of the Geofit input files
const GEOFITFILES =[
  "in_alt.dat",
  "in_lat.dat",
  "in_press.dat",
  "in_temp.dat",
  "in_vrm_prof.dat",
]

# headers of the Geofit input files
const GEOFITHEADERS =[
  24, # in_alt.dat
  21, # in_lat.dat
  19, # in_press.dat
  19, # in_temp.dat
  20, # in_vrm_prof.dat
]
