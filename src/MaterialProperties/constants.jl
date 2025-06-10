using StaticArrays
const CURRENTMATERIALMODEL = ("Cidor","Mathar")
############################################################################################
# CONSTANTS FOR CONVERSION
############################################################################################# TEMPERATURE CONVERSIONS
const CONVERSIONKELVINTOCELSIUS::Float64 = convert(Float64,ustrip(uconvert(°C,0*K))) # °C
const CONVERSIONCELSIUSTOKELVIN::Float64 = -CONVERSIONKELVINTOCELSIUS # K
# PRESSURE CONVERSIONS
const CONVERSIONPASCALTOHPA::Float64     = convert(Float64,ustrip(uconvert(hPa,1*Pa))) # hPa
const CONVERSIONPASCALTOMBAR::Float64    = CONVERSIONPASCALTOHPA # mba
const CONVERSIONPASCALTOATM::Float64    = convert(Float64,ustrip(uconvert(atm,1*Pa))) # atm

const CONVERSIONHPATOPASCAL::Float64     = 1/CONVERSIONPASCALTOHPA # Pa
const CONVERSIONMBARTOPASCAL::Float64    = 1/CONVERSIONPASCALTOMBAR # Pa
const CONVERSIONATMTOMBAR::Float64       = 1/CONVERSIONPASCALTOATM # Pa

const CONVERSIONCMTOμM::Float64 = convert(Float64,ustrip(uconvert(u"μm",1*u"cm"))) # cm
const CONVERSIONμMTOCM::Float64 = 1/CONVERSIONCMTOμM # μm
###########################################################################################
# CIDDOR MODEL CONSTANTS 1996
#############################################################################################
const CID_K0::Float64 = convert(Float64,238.0185) # µm^2
const CID_K1::Float64 = convert(Float64,5792105)  # µm^2
const CID_K2::Float64 = convert(Float64,57.362)   # µm^2
const CID_K3::Float64 = convert(Float64,167917)   # µm^2
const CIDREF_CO2::Float64 = convert(Float64,0.534e-6) # unitless
const CIDREF_CO2PPM::Float64 = convert(Float64,450) # ppm
const CID_WVCORRECTING::Float64 = convert(Float64,1.022e-8) # unitless
const CID_W0::Float64 = convert(Float64,295.235)    # μm^-2 in the paper but it should be unitless
const CID_W1::Float64 = convert(Float64,2.6422)     # μm^2
const CID_W2::Float64 = convert(Float64,-0.032380)   # μm^4
const CID_W3::Float64 = convert(Float64,0.004028)  # μm^6

const CID_A::Float64 = convert(Float64,1.2378847e-5)  #K^-2
const CID_B::Float64 = convert(Float64,-1.9121316e-2) #K^-1
const CID_C::Float64 = convert(Float64,33.93711047)
const CID_D::Float64 = convert(Float64,-6.3431645e3)  #K
const CID_E::Float64 = convert(Float64,-2663.5) #K
const CID_F::Float64 =  12.537

const CID_α::Float64 = convert(Float64,1.00062)
const CID_β::Float64 = convert(Float64,3.14e-8)       #Pa^-1,
const CID_γ::Float64 = convert(Float64,5.6e-7)        #°C^-2
const CIDREF_MOLARMASSW = Mw::Float64 = convert(Float64,0.018015)
const CIDREF_MOLARMASSAIR::Float64 = convert(Float64,28.9635) # g/mol
const CIDREF_MOLARMASSCO2::Float64 = convert(Float64,12.011e-6) # g/mol
const CIDREF_MOLARMASSC02PPM::Float64 = convert(Float64,400) # ppm
const RGAS::Float64 = convert(Float64,8.314510) # J/(mol K) gas constant
const CIDREF_PRESSUREAIR::Float64 = convert(Float64,101325) # Pa
const CIDREF_PRESSUREWATER::Float64 = convert(Float64,1333) # Pa
const CIDREF_TEMPERATUREKELVINAIR::Float64 = convert(Float64,288.15) # K
const CIDREF_TEMPERATURECELSIUSAIR::Float64 = ustrip(uconvert(°C,CIDREF_TEMPERATUREKELVINAIR*K)) # °C
const CIDREF_TEMPERATUREWATERKELVIN::Float64 = convert(Float64,293.15) # K
const CIDREF_TEMPERATUREWATERCELSIUS::Float64 = ustrip(uconvert(°C,CIDREF_TEMPERATUREWATERKELVIN*K)) # °C

const CID_COMPRESSIBILITYA0::Float64 = convert(Float64,1.58123e-6) # K/Pa
const CID_COMPRESSIBILITYA1::Float64 = convert(Float64,-2.9331e-8) # 1/Pa
const CID_COMPRESSIBILITYA2::Float64 = convert(Float64,1.1043e-10) # 1/K/Pa
const CID_COMPRESSIBILITYB0::Float64 = convert(Float64,5.707e-6) # K/Pa
const CID_COMPRESSIBILITYB1::Float64 = convert(Float64,-2.051e-8) # 1/Pa
const CID_COMPRESSIBILITYC0::Float64 = convert(Float64,1.9898e-4) # K/Pa
const CID_COMPRESSIBILITYC1::Float64 = convert(Float64,-2.376e-6) # 1/Pa
const CID_COMPRESSIBILITYD::Float64 = convert(Float64,1.83e-11) # K^2/Pa^-2
const CID_COMPRESSIBILITYE::Float64 = convert(Float64,-0.765e-8) # K^2/Pa^-2
#########################################################################################################################

#########################################################################################################################
# CARLOTTI MODEL CONSTANTS
#########################################################################################################################
const CARLOTTI_BASE=0.000272632                      # unitless, base constant for the Carlotti model
const CARLOTTIREF_PRESSTEMP=288.16/101325       # K/Pa if multiplied by P/T standart 15°C and 101325 Pa returns 1.0
#########################################################################################################################


#########################################################################################################################
# MATHAR MODEL CONSTANTS
#########################################################################################################################
const MATREF_T = 273.15+17.5 # K
const MATREF_P = 75000       # Pa
const MATREF_H = 10.0 # % relative humidity
#########################################################################################################################
# M013_025  valid in the range of 1.3μm - 2.5μm
#########################################################################################################################cref = [ 0.200192e-3,  0.113474e-9,  -0.424595e-14,  0.100957e-16, -0.293315e-20,  0.307228e-24] # cm^j
const M013_025REFC = [ 0.200192e-3,  0.113474e-9,  -0.424595e-14,  0.100957e-16, -0.293315e-20,  0.307228e-24] # cm^j

const M013_025CT   = [ 0.588625e-1, -0.385766e-7,   0.888019e-10, -0.567650e-13,  0.166615e-16, -0.174845e-20] # cm^j · K
const M013_025CTT  = [-3.01513,      0.406167e-3,  -0.514544e-6,   0.343161e-9,  -0.101189e-12,  0.106749e-16] # cm^j · K^2
const M013_025CH   = [-0.103945e-7,  0.136858e-11, -0.171039e-14,  0.112908e-17, -0.329925e-21,  0.344747e-25] # cm^j · %^-1
const M013_025CHH  = [ 0.573256e-12, 0.186367e-16, -0.228150e-19,  0.150947e-22, -0.441214e-26,  0.461209e-30] # cm^j · %^-2
const M013_025CP   = [ 0.267085e-8,  0.135941e-14,  0.135295e-18,  0.818218e-23, -0.222957e-26,  0.249964e-30] # cm^j · Pa^-1
const M013_025CPP  = [ 0.609186e-17, 0.519024e-23, -0.419477e-27,  0.434120e-30, -0.122445e-33,  0.134816e-37] # cm^j · Pa^-2
const M013_025CTH  = [ 0.497859e-4, -0.661752e-8,   0.832034e-11, -0.551793e-14,  0.161899e-17, -0.169901e-21] # cm^j · K · %^-1
const M013_025CTP  = [ 0.779176e-6,  0.396499e-12,  0.395114e-16,  0.233587e-20, -0.636441e-24,  0.716868e-28] # cm^j · K · Pa^-1
const M013_025CHP  = [-0.206567e-15, 0.106141e-20, -0.149982e-23,  0.984046e-27, -0.288266e-30,  0.299105e-34] # cm^j · %^-1 · Pa^-1
# reference wavelength,
const M013_025REFσ = CONVERSIONCMTOμM/2.25    # cm^−1
#########################################################################################################################
# M028_042  valid in the range of 2.8μm - 4.2μm
#########################################################################################################################const M028_042REFC = @SVector [ 0.199885e-3,  0.344739e-9,  -0.273714e-12,  0.393383e-15, -0.569488e-17,  0.164556e-19] # cm^j
const M028_042REFC  = @SVector [ 0.200049e-3,  0.145221e-9,   0.250951e-12, -0.745834e-15, -0.161432e-17,  0.352780e-20] # cm^j
const M028_042CT    = @SVector  [ 0.588432e-1, -0.825182e-7,   0.137982e-9,   0.352420e-13, -0.730651e-15, -0.167911e-18] # cm^j · K
const M028_042CTT   = @SVector  [-3.13579,      0.694124e-3,  -0.500604e-6,  -0.116668e-8,   0.209644e-11,  0.591037e-14] # cm^j · K^2
const M028_042CH    = @SVector  [-0.108142e-7,  0.230102e-11, -0.154652e-14, -0.323014e-17,  0.630616e-20,  0.173880e-22] # cm^j · %^-1
const M028_042CHH   =  @SVector [ 0.586812e-12, 0.312198e-16, -0.197792e-19, -0.461945e-22,  0.788398e-25,  0.245580e-27] # cm^j · %^-2
const M028_042CP    =  @SVector [ 0.266900e-8,  0.168162e-14,  0.353075e-17, -0.963455e-20, -0.223079e-22,  0.453166e-25] # cm^j · Pa^-1
const M028_042CPP   =  @SVector [ 0.608860e-17, 0.461560e-22,  0.184282e-24, -0.524471e-27, -0.121299e-29,  0.246512e-32] # cm^j · Pa^-2
const M028_042CTH   =  @SVector [ 0.517962e-4, -0.112149e-7,   0.776507e-11,  0.172569e-13, -0.320582e-16, -0.899435e-19] # cm^j · K · %^-1
const M028_042CTP   =  @SVector [ 0.778638e-6,  0.446396e-12,  0.784600e-15, -0.195151e-17, -0.542083e-20,  0.103530e-22] # cm^j · K · Pa^-1
const M028_042CHP   = @SVector [-0.217243e-15, 0.104747e-20, -0.523689e-23,  0.817386e-26,  0.309913e-28, -0.363491e-31] # cm^j · %^-1 · Pa^-1
# reference wavelength,
const M028_042REFσ       = CONVERSIONCMTOμM/3.4    # cm^−1
#########################################################################################################################
# M043_052  valid in the range of 4.3μm - 5.2μm
#########################################################################################################################
const M043_052REFC = [ 0.200020e-3,  0.275346e-9,   0.325702e-12, -0.693603e-14,  0.285610e-17,  0.338758e-18] # cm^j
const M043_052CT   = [ 0.590035e-1, -0.375764e-6,   0.134585e-9,   0.124316e-11,  0.508510e-13, -0.189245e-15] # cm^j · K
const M043_052CTT  = [-4.09830,      0.250037e-2,   0.275187e-6,  -0.653398e-8,  -0.310589e-9,   0.127747e-11] # cm^j · K^2
const M043_052CH   = [-0.140463e-7,  0.839350e-11, -0.190929e-14, -0.121399e-16, -0.898863e-18,  0.364662e-20] # cm^j · %^-1
const M043_052CHH  = [ 0.543605e-12, 0.112802e-15, -0.229979e-19, -0.191450e-21, -0.120352e-22,  0.500955e-25] # cm^j · %^-2
const M043_052CP   = [ 0.266898e-8,  0.273629e-14,  0.463466e-17, -0.916894e-23,  0.136685e-21,  0.413687e-23] # cm^j · Pa^-1
const M043_052CPP  = [ 0.610706e-17, 0.116620e-21,  0.244736e-24, -0.497682e-26,  0.742024e-29,  0.224625e-30] # cm^j · Pa^-2
const M043_052CTH  = [ 0.674488e-4, -0.406775e-7,   0.289063e-11,  0.819898e-13,  0.468386e-14, -0.191182e-16] # cm^j · K · %^-1
const M043_052CTP  = [ 0.778627e-6,  0.593296e-12,  0.145042e-14,  0.489815e-17,  0.327941e-19,  0.128020e-21] # cm^j · K · Pa^-1
const M043_052CHP  = [-0.211676e-15, 0.487921e-20, -0.682545e-23,  0.942802e-25, -0.946422e-27, -0.153682e-29] # cm^j · %^-1 · Pa^-1

# reference wavelength
const M043_052REFσ = CONVERSIONCMTOμM/4.8    # cm^−1
#########################################################################################################################
# M075_141  valid in the range of 7.5μm - 14.1μm  % independent of CO2 concentration
#########################################################################################################################
const M75_141REFC = @SVector [ 0.199885e-3,  0.344739e-9,  -0.273714e-12,  0.393383e-15, -0.569488e-17,  0.164556e-19] # cm^j
const M075_141CT         = @SVector [ 0.593900e-1, -0.172226e-5,   0.237654e-8,  -0.381812e-11,  0.305050e-14, -0.157464e-16] # cm^j · K
const M075_141CTT        = @SVector [-6.50355,      0.103830e-1,  -0.139464e-4,   0.220077e-7,  -0.272412e-10,  0.126364e-12] # cm^j · K^2
const M075_141CH         = @SVector  [-0.221938e-7,  0.347377e-10, -0.465991e-13,  0.735848e-16, -0.897119e-19,  0.380817e-21] # cm^j · %^-1
const M075_141CHH        = @SVector  [ 0.393524e-12, 0.464083e-15, -0.621764e-18,  0.981126e-21, -0.121384e-23,  0.515111e-26] # cm^j · %^-2
const M075_141CP         =  @SVector [ 0.266809e-8,  0.695247e-15,  0.159070e-17, -0.303451e-20, -0.661489e-22,  0.178226e-24] # cm^j · Pa^-1
const M075_141CPP        =  @SVector [ 0.610508e-17, 0.227694e-22,  0.786323e-25, -0.174448e-27, -0.359791e-29,  0.978307e-32] # cm^j · Pa^-2
const M075_141CTH        =  @SVector [ 0.106776e-3, -0.168516e-6,   0.226201e-9,  -0.356457e-12,  0.437980e-15, -0.194545e-17] # cm^j · K · %^-1
const M075_141CTP        =  @SVector [ 0.77368e-6,   0.216404e-12,  0.581805e-15, -0.189618e-17, -0.198869e-19,  0.589381e-22] # cm^j · K · Pa^-1
const M075_141CHP        = @SVector  [-0.206365e-15, 0.300234e-19, -0.426519e-22,  0.684306e-25, -0.467320e-29,  0.126117e-30] # cm^j · %^-1 · Pa^-1
# reference wavelength,
const M075_141REFσ       = CONVERSIONCMTOμM/10.1    # cm^−1
#########################################################################################################################
# M160_240  valid in the range of 16.0μm - 24.0μm
#########################################################################################################################
const M160_240REFC = [ 0.199436e-3,  0.299123e-8,  -0.214862e-10,  0.143338e-12,  0.122398e-14, -0.114628e-16] # cm^j
# something seems to be wrong with cT...
const M160_240CT   = [ 0.621723e-1, -0.177074e-4,   0.152213e-6,  -0.954584-9,   -0.996706e-11,  0.921476e-13] # cm^j · K
const M160_240CTT  = [-23.2409,      0.108557,     -0.102439e-2,   0.634072e-5,   0.762517e-7,  -0.675587e-9 ] # cm^j · K^2
const M160_240CH   = [-0.772707e-7,  0.347237e-9,  -0.272675e-11,  0.170858e-13,  0.156889e-15, -0.150004e-17] # cm^j · %^-1
const M160_240CHH  = [-0.326604e-12, 0.463606e-14, -0.364272e-16,  0.228756e-18,  0.209502e-20, -0.200547e-22] # cm^j · %^-2
const M160_240CP   = [ 0.266827e-8,  0.120788e-14,  0.522646e-17,  0.783027e-19,  0.753235e-21, -0.228819e-24] # cm^j · Pa^-1
const M160_240CPP  = [ 0.613675e-17, 0.585494e-22,  0.286055e-24,  0.425193e-26,  0.413455e-28, -0.812941e-32] # cm^j · Pa^-2
const M160_240CTH  = [ 0.375974e-3, -0.171849e-5,   0.146704e-7,  -0.917231e-10, -0.955922e-12,  0.880502e-14] # cm^j · K · %^-1
const M160_240CTP  = [ 0.778436e-6,  0.461840e-12,  0.306229e-14, -0.623183e-16, -0.161119e-18,  0.800756e-20] # cm^j · K · Pa^-1
const M160_240CHP  = [-0.272614e-15, 0.304662e-18, -0.239590e-20,  0.149285e-22,  0.136086e-24, -0.130999e-26] # cm^j · %^-1 · Pa^-1
# reference wavelength,
const M160_240REFσ = CONVERSIONCMTOμM/20.0    # cm^−1
