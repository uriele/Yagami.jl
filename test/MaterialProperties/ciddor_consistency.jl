using Test
using Yagami.MaterialProperties
using Yagami.MaterialProperties:celsius_to_kelvin
using Yagami.MaterialProperties:CIDREF_PRESSUREAIR
using Yagami
using PyCall
pyimport("numpy")

# note: I use testset to keep the scope of the PyCall variables local
@testset "Ciddor consistency" begin

py"""
# Crated by Mikhail Polyanskiy https://github.com/polyanskiy
# This script is adapted from the original
# Original data: Ciddor 1996, https://doi.org/10.1364/AO.35.001566
# CHANGELOG
# 2017-11-23 [Misha Polyanskiy] original version
# 2024-03-03 [Misha Polyanskiy] minor refactoringpy"""

import numpy as np
π = np.pi


def Z(T,p,xw):
    t=T-273.15
    a0 = 1.58123e-6   #K·Pa^-1
    a1 = -2.9331e-8   #Pa^-1
    a2 = 1.1043e-10   #K^-1·Pa^-1
    b0 = 5.707e-6     #K·Pa^-1
    b1 = -2.051e-8    #Pa^-1
    c0 = 1.9898e-4    #K·Pa^-1
    c1 = -2.376e-6    #Pa^-1
    d  = 1.83e-11     #K^2·Pa^-2
    e  = -0.765e-8    #K^2·Pa^-2
    return 1-(p/T)*(a0+a1*t+a2*t**2+(b0+b1*t)*xw+(c0+c1*t)*xw**2) + (p/T)**2*(d+e*xw**2)


def n(λ,t,p,h,xc):
    # λ: wavelength, 0.3 to 1.69 μm
    # t: temperature, -40 to +100 °C
    # p: pressure, 80000 to 120000 Pa
    # h: fractional humidity, 0 to 1
    # xc: CO2 concentration, 0 to 2000 ppm

    σ = 1/λ           #μm^-1

    T= t + 273.15     #Temperature °C -> K

    R = 8.314510      #gas constant, J/(mol·K)

    k0 = 238.0185     #μm^-2
    k1 = 5792105      #μm^-2
    k2 = 57.362       #μm^-2
    k3 = 167917       #μm^-2

    w0 = 295.235      #μm^-2
    w1 = 2.6422       #μm^-2
    w2 = -0.032380    #μm^-4
    w3 = 0.004028     #μm^-6

    A = 1.2378847e-5  #K^-2
    B = -1.9121316e-2 #K^-1
    C = 33.93711047
    D = -6.3431645e3  #K

    α = 1.00062
    β = 3.14e-8       #Pa^-1,
    γ = 5.6e-7        #°C^-2

    #saturation vapor pressure of water vapor in air at temperature T (Pa)
    svp = np.where(t>=0,
        np.exp(A*T**2 + B*T + C + D/T), # if t>=0
        10**(-2663.5/T+12.537))         # if t<0

    #enhancement factor of water vapor in air
    f = α + β*p + γ*t**2

    #molar fraction of water vapor in moist air
    xw = f*h*svp/p

    #refractive index of standard air at 15 °C, 101325 Pa, 0% humidity, 450 ppm CO2
    nas = 1 + (k1/(k0-σ**2)+k3/(k2-σ**2))*1e-8

    #refractive index of standard air at 15 °C, 101325 Pa, 0% humidity, xc ppm CO2
    naxs = 1 + (nas-1) * (1+0.534e-6*(xc-450))

    #refractive index of water vapor at standard conditions (20 °C, 1333 Pa)
    nws = 1 + 1.022*(w0+w1*σ**2+w2*σ**4+w3*σ**6)*1e-8

    Ma = 1e-3*(28.9635 + 12.011e-6*(xc-400)) #molar mass of dry air, kg/mol
    Mw = 0.018015                            #molar mass of water vapor, kg/mol

    Za = Z(288.15, 101325, 0)                #compressibility of dry air
    Zw = Z(293.15, 1333, 1)                  #compressibility of pure water vapor

    #Eq.4 with (T,P,xw) = (288.15, 101325, 0)
    ρaxs = 101325*Ma/(Za*R*288.15)           #density of standard air

    #Eq 4 with (T,P,xw) = (293.15, 1333, 1)
    ρws  = 1333*Mw/(Zw*R*293.15)             #density of standard water vapor

    # two parts of Eq.4: ρ=ρa+ρw
    ρa   = p*Ma/(Z(T,p,xw)*R*T)*(1-xw)       #density of the dry component of the moist air
    ρw   = p*Mw/(Z(T,p,xw)*R*T)*xw           #density of the water vapor component

    nprop = 1 + (ρa/ρaxs)*(naxs-1) + (ρw/ρws)*(nws-1)

    return nprop

def ff(t,p):
  α = 1.00062
  β = 3.14e-8       #Pa^-1,
  γ = 5.6e-7        #°C^-2
  α + β*p + γ*t**2

def xm(t,p,h):
    f= ff(t,p)
    svp = sanity_check1(t)
    xw = f*h*svp/p

def sanity_check1(t):
    A = 1.2378847e-5  #K^-2
    B = -1.9121316e-2 #K^-1
    C = 33.93711047
    D = -6.3431645e3  #K
    T= t + 273.15     #Temperature °C -> K
    return np.where(t>=0,
        np.exp(A*T**2 + B*T + C + D/T), # if t>=0
        10**(-2663.5/T+12.537))
"""

n_py = py"n"
sanity_check_py1 = py"sanity_check1"
xm_py = py"xm"
ff_py = py"ff"

λ= LinRange(0.3,1.69,10) # wavelength in µm
λ=vcat(λ...,LinRange(9.5,11.0,10)...)

t= LinRange(-40,100,50) # temperature in °C
T= @. celsius_to_kelvin(t) # temperature in K
p= LinRange(800,CIDREF_PRESSUREAIR,50) # pressure in Pa
h= LinRange(0,1,10) # humidity in fraction
co2ppm= LinRange(0,500,20) # CO2 concentration in ppm
n_mine = Array{Float64}(undef,length(λ),length(t),length(p),length(h),length(co2ppm));
n_python = similar(n_mine);
svp_mine = similar(n_mine);
svp_python = similar(n_mine);

# correct for really small values since the equation is not continuous across 0
t=collect(t)
t[abs.(t).< 1e-10].=0.0
################################################################################

@test(__saturation_vapor_pressure.(T)≈ sanity_check_py1(t))


for (i,(t,T)) in enumerate(zip(t,T))
   for (j,p) in enumerate(p)
     for (k,h) in enumerate(h)
      for (l,co2ppm) in enumerate(co2ppm)
          @inbounds n_mine[:,i,j,k,l]      .= @. Yagami.MaterialProperties.ciddor_refractive_index(T,p,λ,h,co2ppm)
          @inbounds n_python[:,i,j,k,l]    .= n_py(λ,t,p,h,co2ppm)
          @inbounds svp_mine[:,i,j,k,l]    .=@. Yagami.MaterialProperties.__saturation_vapor_pressure(T)
          @inbounds svp_python[:,i,j,k,l]  .= sanity_check_py1(t)
      end
    end
  end
end

@test svp_mine ≈ svp_python
@test n_mine ≈ n_python
