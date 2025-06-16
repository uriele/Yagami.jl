
###############################################################################
# Constants for the numerical methods
###############################################################################
# Complement of the golden ratio
const CGOLD = 1-1/Base.MathConstants.golden
const GLIMIT = 100 # maximum allowed step size in the parabolic fit
const GOLDEN = Base.MathConstants.golden
const RGOLD = Base.MathConstants.golden-1
const TOL= 1e-10
############################################################################################
# STANDARD EARTH ELLIPSOID
############################################################################################
const WGS84EARTH    = ellipsoid(WGS84Latest)
###
const WGS84MAJORAXIS = majoraxis(WGS84EARTH) |> er-> uconvert(km,er) |> ustrip
const WGS84MINORAXIS = minoraxis(WGS84EARTH) |> er-> uconvert(km,er) |> ustrip
const WGS84ECCENTRICITY² = eccentricity²(WGS84EARTH) |> ustrip
const WGS84ECCENTRICITY = sqrt(WGS84ECCENTRICITY²)   |> ustrip
const WGS84NORMALIZEMINORAXIS = WGS84MINORAXIS / WGS84MAJORAXIS
const WGS84COMPLECCENTRICITY² = 1 - WGS84ECCENTRICITY²
###
const REFMAJORAXIS = Ref{Float64}(WGS84MAJORAXIS)
const REFMINORAXIS = Ref{Float64}(WGS84MINORAXIS)
const REFECCENTRICITY² = Ref{Float64}(WGS84ECCENTRICITY²)
const REFECCENTRICITY = Ref{Float64}(WGS84ECCENTRICITY)
const REFNORMALIZEMINORAXIS = Ref{Float64}(WGS84NORMALIZEMINORAXIS)
const REFCOMPLECCENTRICITY² = Ref{Float64}(WGS84COMPLECCENTRICITY²)

############################################################################################
const DATUMINFO = (:MAJORAXIS, :MINORAXIS, :ECCENTRICITY², :ECCENTRICITY, :NORMALIZEMINORAXIS, :COMPLECCENTRICITY²)
const SETDATUM  = Expr[
    :(REFMAJORAXIS[]          = _MAJORAXIS),
    :(REFMINORAXIS[]          = _MINORAXIS),
    :(REFECCENTRICITY²[]      = _ECCENTRICITY²),
    :(REFECCENTRICITY[]       = _ECCENTRICITY),
    :(REFNORMALIZEMINORAXIS[] = _NORMALIZEMINORAXIS),
    :(REFCOMPLECCENTRICITY²[] = _COMPLECCENTRICITY²)
]
############################################################################################
