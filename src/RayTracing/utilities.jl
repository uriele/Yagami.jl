gettemperature(atm::AtmosphereSetting, θ, h) = atm.temperature(θ, h)
getpressure(atm::AtmosphereSetting, θ, h) = atm.pressure(θ, h)
gethumidity(atm::AtmosphereSetting, θ, h) = atm.humidity(θ, h)
getco2ppm(atm::AtmosphereSetting, θ, h) = atm.co2ppm(θ, h)
getwavelength(atm::AtmosphereSetting, θ, h) = atm.wavelength(θ, h)

getknotsh(atm::AtmosphereSetting) = atm.knots_h
getknotsθ(atm::AtmosphereSetting) = atm.knots_θ
