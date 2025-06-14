
# Type of interpolation used in atmospheric models.
abstract type AtmosphericInterpolation end
struct BiLinear<:AtmosphericInterpolation  end
struct LogLinear<:AtmosphericInterpolation end

# Type of mean used in atmospheric models.
abstract type MeanType end
struct ArithmeticMean <: MeanType end
struct GeometricMean <: MeanType end
struct LogMean <: MeanType end


const OFFSET0 = 1e-10



"""
  `AtmInteoplate{N,M,T}`

A structure representing an atmospheric interpolation object for ray tracing.

  `AtmInterpolate(knots_θ, knots_h, parameters; logh=false)`

Create an atmospheric interpolation object for ray tracing. The interpolation is defined by the knots of height and angle, and the parameters matrix.
  # Arguments
  - `knots_θ::AbstractVector{T}`: A vector of angles in degrees, where `T` is a subtype of `Real`.
  - `knots_h::AbstractVector{T}`: A vector of heights in kilometers, where `T` is a subtype of `Real`.
  - `parameters::AbstractMatrix{T}`: A matrix of parameters for the interpolation, where `T` is a subtype of `Real`.
  - `logh::Bool`: If true, use logarithmic interpolation for height. Default is false.
  # Returns
  - An `AtmInterpolate` object with the specified knots and parameters.
"""
struct AtmInterpolate{N,M,T<:AbstractFloat,VN<:AbstractVector{T},VM<:AbstractVector{T},MP<:AbstractMatrix{T},ITP<:AtmosphericInterpolation}
  knots_h::VM
  knots_θ::VN
  A ::MP
  __θmin::T
  __θmax::T
  __hmin::T
  __hmax::T


  function AtmInterpolate(knots_θ::AbstractVector{T},knots_h::AbstractVector{T},
     parameters::AbstractMatrix{T};logh::Bool=false) where {T}

    @inline isunique(x::AbstractVector) = length(x) == length(unique(x))

    # initialize elements of AtmInterpolate
    knots_θ = collect(knots_θ)
    knots_h = collect(knots_h)
    ########################################
    ITP = logh ? LogLinear : BiLinear
    N,M= size(parameters)
    Ni,Mi=length(knots_θ),length(knots_h)

    @assert(N == Ni && M == Mi, "knots_θ and knots_h must have the same length as the number of rows and columns of parameters but got knots_θ=$knots_θ, knots_h=$knots_h, parameters=$parameters.")
    ##############################################################################################################
    @assert(issorted(knots_h, rev=true), "knots_h must be sorted in descending order but got knots_h=$knots_h.")
    @assert(minimum(knots_h) >= 0, "knots_h must be non-negative but got knots_h=$knots_h.")
    @assert(isunique(knots_h), "knots_h must be unique but got knots_h=$knots_h.")
    ###############################################################################################################
    knots_θ .= @. mod(knots_θ, 360)
    @assert(isunique(knots_θ), "knots_θ must be unique but got knots_θ=$knots_θ.")
    @assert(issorted(knots_θ), "knots_θ must be sorted in ascending order in the range [0,360).") # this would help with generating the
    ###############################################################################################################
    # interpolation matrix
    A=vcat(parameters,parameters[1:1,:]) # add the first column to the end to make it periodic
    __θmin = knots_θ[1]
    __θmax = knots_θ[end]
    __hmin = knots_h[end]
    __hmax = knots_h[1]

    knots_θ = vcat(knots_θ,knots_θ[1]) # remove the last element since it is the same as the first one

    VN=typeof(knots_θ)
    VM=typeof(knots_h)
    MP=typeof(A)
    return new{N,M,T,VN,VM,MP,ITP}(knots_h, knots_θ, A,__θmin, __θmax, __hmin, __hmax)
  end
end

Base.propertynames(::AtmInterpolate) = (:knots_h, :knots_θ, :N,:M,:interpolation)


function Base.show(io::IO,atm::AtmInterpolate)
  __θmin= round(getfield(atm,:__θmin), digits=3)
  __θmax= round(getfield(atm,:__θmax), digits=3)
  __hmin= round(getfield(atm,:__hmin), digits=3)
  __hmax= round(getfield(atm,:__hmax), digits=3)
  print(io, "Atm([",__θmin,"°..",__θmax,"°]], [",__hmin,"km..",__hmax,"km]")
end

function Base.show(io::IO,::Type{<:AtmInterpolate{N,M,T,VN,VM,MP,ITP}}) where {N,M,T<:AbstractFloat,VN,VM,MP,ITP}
  print(io, "AtmInterpolate{",N,",",M,",",T,"} (with $(ITP==BiLinear ? "BiLinear" : "LogLinear") interpolation)")
end

function Base.show(io::IO, ::MIME"text/plain",atm::AtmInterpolate{N,M,T,VN,VM,MP,ITP}) where {N,M,T<:AbstractFloat,VN,VM,MP,ITP}
  __θmin= round(getfield(atm,:__θmin), digits=3)
  __θmax= round(getfield(atm,:__θmax), digits=3)
  __hmin= round(getfield(atm,:__hmin), digits=3)
  __hmax= round(getfield(atm,:__hmax), digits=3)

  println(io, "AtmInterpolate{",N,",",M,",", T,"}")
  println(io, "  knots_θ: $(__θmin)°..$(__θmax)°")
  println(io, "  knots_h: $(__hmin)km..$(__hmax)km")
  println(io, "  Interpolation: $(ITP == BiLinear ? "BiLinear" : "LogLinear")")
end


function Base.getproperty(atm::AtmInterpolate{N,M,T,VN,VM,MP,ITP}, sym::Symbol) where {N,M,T<:AbstractFloat,VN,VM,MP,ITP}
  if sym == :knots_h
    return getfield(atm,:knots_h)
  elseif sym == :knots_θ
    return getfield(atm,:knots_θ)
  elseif sym == :interpolation
    return ITP
  elseif sym == :N
    return N
  elseif sym == :M
    return M
  end
end

function Base.setproperty!(::AtmInterpolate, ::Symbol, value)
  error("Cannot reset the function of the BracketFuncWrapper after creation.")
end



@inline __geth(h::T,knots_h::AbstractVector{T},__hmin::T,__hmax::T,M::Int) where {T<:Real} = if __hmin< h < __hmax
  searchsortedlast(knots_h, h,rev=true)
  elseif h == knots_h[1]
    1
  elseif h == knots_h[end]
    M
  else
    throw(BoundsError("h must be in the range [$__hmin, $__hmax]."))
  end

@inline __getθ(θ::T,knots_θ::AbstractVector{T},__θmin::T,__θmax::T,N::Int) where {T<:Real} = begin
  θ = mod(θ, 360);
  if __θmin <= θ < __θmax
    searchsortedlast(knots_θ, θ)
  else
    N
  end
end


@inline (itp::AtmInterpolate{N,M,T,VN,VM,MP,ITP})(θ,h) where {N,M,T<:Real,VM<:AbstractVector{T},VN<:AbstractVector{T},MP<:AbstractMatrix{T},ITP<:BiLinear} = begin

  __θmin = getfield(itp, :__θmin)
  __θmax = getfield(itp, :__θmax)
  __hmin = getfield(itp, :__hmin)
  __hmax = getfield(itp, :__hmax)

  θ=mod(θ, 360) # ensure θ is in the range [0,360)
  j = __geth(h, itp.knots_h, __hmin, __hmax, M)
  i = __getθ(θ, itp.knots_θ, __θmin, __θmax, N)

  Δθ = mod(itp.knots_θ[i+1] - itp.knots_θ[i],360)
  fθ = mod(θ - itp.knots_θ[i],360) / Δθ



  aij=(getfield(itp,:A)[i,j]*(1-fθ)+(fθ)*getfield(itp,:A)[i+1,j])


  if j==M
    return aij
  end


  aij1=(getfield(itp,:A)[i,j+1]*(1-fθ)+(fθ)*getfield(itp,:A)[i+1,j+1])

  Δh = (itp.knots_h[j+1] - itp.knots_h[j])
  fh = (h - itp.knots_h[j]) / Δh



  return aij1*(fh)+aij*(1-fh)
end


@inline (itp::AtmInterpolate{N,M,T,VM,VN,MP,ITP})(θ,h) where {N,M,T<:Real,VM<:AbstractVector{T},VN<:AbstractVector{T},MP<:AbstractMatrix{T},ITP<:LogLinear} = begin

  __θmin = getfield(itp, :__θmin)
  __θmax = getfield(itp, :__θmax)
  __hmin = getfield(itp, :__hmin)
  __hmax = getfield(itp, :__hmax)
  θ=mod(θ, 360) # ensure θ is in the range [0,360)
  j = __geth(h, itp.knots_h, __hmin, __hmax, M)
  i = __getθ(θ, itp.knots_θ, __θmin, __θmax, N)

  Δθ = mod(itp.knots_θ[i+1] - itp.knots_θ[i],360)
  fθ = mod(θ - itp.knots_θ[i],360) / Δθ



  aij=(getfield(itp,:A)[i,j]*(1-fθ)+(fθ)*getfield(itp,:A)[i+1,j])


  if j==M
    return aij
  end


  aij1=(getfield(itp,:A)[i,j+1]*(1-fθ)+(fθ)*getfield(itp,:A)[i+1,j+1])

  Δh = (itp.knots_h[j+1] - itp.knots_h[j])
  fh = (h - itp.knots_h[j]) / Δh



  return aij1^(fh)*aij^(1-fh)
end




@inline __mean(::MeanType,topleft::T,topright::T,bottomleft::T,bottomright::T)  where T = error("Mean type not implemented.")
@inline __mean(::ArithmeticMean, topleft::T, topright::T, bottomleft::T, bottomright::T) where T = (topleft + topright + bottomleft + bottomright) / 4
@inline function __mean(::GeometricMean, topleft::T, topright::T, bottomleft::T, bottomright::T) where T
  top    =(topleft+topright)/2
  bottom =(bottomleft+bottomright)/2
  exp(log(top*bottom)/2)
end

@inline function __mean(::LogMean, topleft::T, topright::T, bottomleft::T, bottomright::T) where T
  top    =(topleft+topright)/2
  bottom =(bottomleft+bottomright)/2
  top==bottom ? top : (top-bottom)/(log(top/bottom))
end
