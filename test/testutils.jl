

_N₀(θ,majoraxis=WGS84MAJORAXIS,minoraxis=WGS84MINORAXIS) = majoraxis^2/sqrt(majoraxis^2*cosd(θ)^2+minoraxis^2*sind(θ)^2)

# Compute the ellipsiod XY exacly at angle θ and height
ellipsepoint(h,θ,majoraxis=WGS84MAJORAXIS,minoraxis=WGS84MINORAXIS) = ((_N₀(θ,majoraxis,minoraxis)+h)*cosd(θ), (minoraxis^2/majoraxis^2*_N₀(θ,majoraxis,minoraxis)+h)*sind(θ))


using PyCall
# to load python scripts from file
pyimport("numpy")
pyimport_conda("importlib.util", "importlib")
importlib = pyimport("importlib.util")
function load_python_scripts(filename::String)
  spec = importlib.spec_from_file_location("mymodule", filename)
  mymodule = importlib.module_from_spec(spec)
  spec.loader.exec_module(mymodule)
  return mymodule
end

interpfunc1(θ,h,hmax=120.0,p0=100) = p0*exp(-h/hmax)
interpfunc2(θ,h,hmax=120.0,p0=100) = p0*sind(θ/2)^2
interpfunc3(θ,h,hmax=120.0,p0=100) = interpfunc1(θ,h,hmax,p0) + interpfunc2(θ,h,hmax,p0)
interpfunc4(θ,h,hmax=120.0,p0=100) = interpfunc1(θ,h,hmax,p0) * interpfunc2(θ,h,hmax,p0)

# Syntetic atmosphere functions for testing
temp_func(θ,h) = 20*sind(θ)+ if h<10.0 # troposhpere
  celsius_to_kelvin(20.0)-80/10.0*h
elseif h<=20.0                         # tropopause
  celsius_to_kelvin(-60.)
elseif h<=50                           # stratosphere
  celsius_to_kelvin(-60.0)+60.0/30.0*(h-20.0)
elseif h<= 55                      # stratopause
  celsius_to_kelvin(0.0)
elseif h<=85                           # mesosphere
  celsius_to_kelvin(0.0)-80/(85-55)*(h-55)
elseif h<=100  # termosphere 1
  celsius_to_kelvin(-80.0)+20/(100-85)*(h-85.0)
else # termosphere 2
  celsius_to_kelvin(-60.)+(120)/20*(h-100.0)
end

press_func(θ,h) = atm_to_pascal(1.) * exp(-h/25)
