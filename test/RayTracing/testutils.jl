


_N₀(θ,majoraxis=WGS84MAJORAXIS,minoraxis=WGS84MINORAXIS) = majoraxis^2/sqrt(majoraxis^2*cosd(θ)^2+minoraxis^2*sind(θ)^2)

# Compute the ellipsiod XY exacly at angle θ and height
ellipsepoint(h,θ,majoraxis=WGS84MAJORAXIS,minoraxis=WGS84MINORAXIS) = ((_N₀(θ,majoraxis,minoraxis)+h)*cosd(θ), (minoraxis^2/majoraxis^2*_N₀(θ,majoraxis,minoraxis)+h)*sind(θ))



interpfunc1(θ,h,hmax=120.0,p0=100) = p0*exp(-h/hmax)
interpfunc2(θ,h,hmax=120.0,p0=100) = p0*sind(θ/2)^2
interpfunc3(θ,h,hmax=120.0,p0=100) = interpfunc1(θ,h,hmax,p0) + interpfunc2(θ,h,hmax,p0)
interpfunc4(θ,h,hmax=120.0,p0=100) = interpfunc1(θ,h,hmax,p0) * interpfunc2(θ,h,hmax,p0)
