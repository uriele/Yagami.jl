# Implementation of the Zbrent method for root finding from
# Numerical Recipes 3rd Edition, Section 9.3. pag 454
mutable struct Zbrent{F,T<:AbstractFloat}
    f::F   # function to find the root of
    x::T   # x at the root/minimum
    fx::T  # value of the function at the root/minimum
    a::T   # lower bound of the bracket
    b::T   # upper bound of the bracket
    tol::T # tolerance for convergence
    itermax::Int # maximum number of iterations
    __iter::Int  # number of iterations performed at the last call of findmin
end


"""
    Zbrent{F,T<:AbstractFloat}(f::F, a::T, fx::T, x::T, b::T, tol::T, itermax::Int=100, __iter::Int=0)
Create a new instance of `Zbrent` for finding the minimum of a function `f`.
The parameters are:
- `f`: the function to find the minimum of.
- `itermax`: the maximum number of iterations to perform.
- `tol`: the tolerance for convergence. (note it is not the absolute tolerance)
"""
Zbrent(::Type{<:T},f::F,itermax::Int=100,tol::T=TOL) where {F,T<:AbstractFloat} =
Zbrent{F,T}(f,zero(T),f(zero(T)),zero(T),zero(T),tol,itermax,0)
"""
    `Zbrent([type=Float64]; f::F, itermax::Int=100, tol::T=T(TOL)) where {F,T<:AbstractFloat}`
Create a new instance of `Zbrent` for finding the minimum of a function `f`.
The parameters are:
- `type`: the type of the floating point numbers to use (default is `Float64`).
- `f`: the function to find the minimum of.
- `itermax`: the maximum number of iterations to perform.
- `tol`: the tolerance for convergence (default is `T(TOL)`).
"""
Zbrent(::Type{<:T}=Float64; f::F,itermax::Int=100,tol::T=T(TOL)) where {F,T<:AbstractFloat} =
    Zbrent(T, f, itermax, tol)


# internal helper to set the bracket values
@inline function __setbracket!(z::Zbrent{F,T}, a::T, x::T,b::T) where {F,T<:AbstractFloat}
  setfield!(z, :a, a)
  setfield!(z, :x, x)
  setfield!(z, :fx, z.f(x))
  setfield!(z, :b, b)
end

##################################################################

"""
    findmin(zb::Zbrent{F,T}) where {F,T<:AbstractFloat}
Find the minimum of the function `zb.f` using the Brent's method.
The function `zb.f` must be defined such that it returns a value for a given input `x`.
The method uses a combination of bisection, secant, and inverse quadratic interpolation to find the minimum.
"""
function findmin(zb::Zbrent{F,T}) where {F,T<:AbstractFloat}
    ax = zb.a
    bx = zb.x
    cx = zb.b
    d::T = 0
    e::T = 0
    fmin::T = 0
    ZEPS = eps(T)

    a= (ax<cx ? ax : cx) # lower bound
    b= (ax>cx ? ax : cx) # upper bound

    x=v=w=bx
    fw=fv=fx = zb.fx

    for iter in 0:zb.itermax
      xm =  (a+b)/2                       # midpoint as failsafe method
      tol2 = 2*(tol1=zb.tol*abs(x)+ZEPS)  # set tol1 and tol2
      if abs(x-xm) <= (tol2-(b-a)/2)
        fmin = fx
        zb.x = x
        zb.fx = fx
        zb.__iter = iter
        return
      end

      if abs(e) > tol1 # trial parabolic fit
        r= (x-w)*(fx-fv)
        q= (x-v)*(fx-fw)
        p= (x-v)*q-(x-w)*r;
        q= 2*(q-r)
        if q > 0
          p = -p
        end
        q = abs(q)
        etemp = e
        e = d
        if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
          d=CGOLD*(e=(x >= xm ? a-x : b-x));
        else
          d = p/q
          u = x+d
          if (u-a < tol2 || b-u < tol2)
            d = sign(xm-x)*tol1
          end
        end
      else
        d = CGOLD *(e=(x >= xm ? a-x : b-x))
      end
      u = abs(d) >= tol1 ? x+d : x+sign(d)*tol1  # new point to evaluate
      fu = zb.f(u)  # function value at new point

      if (fu <= fx)  # update minimum
        if (u >= x)
          a = x
        else
          b = x
        end
        v, w, x = w, x, u  # update points
        fv, fw, fx = fw, fx, fu  # update function values

      else  # update bounds
        if (u < x)
          a = u
        else
          b = u
        end

        if (fu <= fw || w == x)
          v,w   = w,u
          fv,fw = fw,fu
        elseif (fu <= fv || v == x || v == w)
          v = u
          fv = fu
        end

      end
    end
    throw(ArgumentError("Maximum number of iterations reached without convergence."))
end

"""
  bracketmin(z::Zbrent{F,T}) where {F,T<:AbstractFloat}
Bracket the function `z.f` to find a minimum. This function sets the initial bounds for the Brent's method.
It uses the golden section search to find a point where the function value is lower than the initial bounds.
"""
function bracketmin(z::Zbrent{F,T}) where {F,T<:AbstractFloat}
  ax=z.tol
  fa = z.f(ax)

  bx= ax+GOLDEN*fa

  fu::T = 0
  fb = z.f(bx)
  # @info "ax: $ax, bx: $bx, fa: $fa, fb: $fb"
  if fb>fa
    ax, bx = bx, ax  # swap if necessary
    fa, fb = fb, fa
  end

  cx=bx+GOLDEN*(bx-ax)  # set cx
  fc = z.f(cx)
  u=cx
  while (fb > fc)  # loop until we find a point with a lower function value
    r = (bx-ax)*(fb-fc)
    q = (bx-cx)*(fb-fa)
    u = bx - ((bx-cx)*q-(bx-ax)*r)/(2*sign(q-r)*max(abs(q-r),eps(T)))
    ulim = bx + GLIMIT*(cx-bx)

    if ((bx-u)*(u-cx) >0)
      fu = z.f(u)  # evaluate function at u
      if (fu < fc)  # update bounds
        ax, bx = bx, u
        fa, fb = fb, fu
        return __setbracket!(z, ax, bx, cx)  # update zbrent bounds
      elseif (fu > fb)  # update bounds
        cx=u
        fc = fu
        # @info "ax: $ax, bx: $bx, cx: $cx, u: $u, fu: $fu"
        return __setbracket!(z, ax, bx, cx)  # update zbrent bounds
      end
      u=cx+ GOLDEN*(cx-bx)  # set u to the next point
      fu = z.f(u)  # evaluate function at u
    elseif ((cx-u)*(u-ulim) > 0)  # u is between cx and ulim
      fu = z.f(u)  # evaluate function at u
      if (fu < fc)  # update bounds
        bx, cx , u = cx, u, u+GOLDEN*(u-cx)  # update ax, bx, and u
        fb, fc , fu = fc, fu, z.f(u)  # update fa, fb, and fu
      end
    elseif ((u-ulim)*(ulim-cx) >= 0)  # limit parabolic u
      u = ulim  # set u to ulim
      fu = z.f(u)  # evaluate function at u
    else  # reject parabolic
      u= cx + GOLDEN*(cx-bx)  # set u to the next point
      fu = z.f(u)  # evaluate function at u
    end
    ax, bx, cx = bx, cx, u  # update ax, bx, and cx
    fa, fb, fc = fb, fc, fu  # update fa, fb, and fc
  end

  __setbracket!(z, ax, bx, cx) # update zbrent bounds


end
