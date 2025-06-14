using Yagami.MaterialProperties: celsius_to_kelvin,atm_to_pascal,kelvin_to_celsius


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
