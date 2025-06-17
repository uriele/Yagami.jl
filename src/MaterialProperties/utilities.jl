
 for from in  ("celsius","kelvin")
   for to  in ("celsius","kelvin")
    if from != to
      fun = Symbol(from, "_to_", to)
      fun! = Symbol(fun, "!")
      CONVERSION = Symbol("CONVERSION", uppercase(string(from)), "TO", uppercase(string(to)))
      @eval $fun(temperature::T) where {T<:AbstractFloat} = temperature + $CONVERSION::T
      @eval function $fun!(temperature::A)::Nothing where {T<:AbstractFloat,A<:AbstractArray{T}}
          @simd for i in eachindex(temperature)
              @inbounds temperature[i] += $CONVERSION
          end
      end
    end
  end
end

for fromto in ("atm","hpa","mbar")
  fun_from = Symbol(fromto, "_to_", :pascal)
  fun_to   = Symbol(:pascal, "_to_", fromto)
  fun_from! = Symbol(fun_from, "!")
  fun_to! = Symbol(fun_to, "!")
  CONVERSIONTOPASCAL = Symbol("CONVERSION", uppercase(fromto), "TOPASCAL")
  CONVERSIONFROMPASCAL = Symbol("CONVERSION", "PASCALTO", uppercase(fromto))
  @eval $fun_from(pressure::T) where T<:AbstractFloat = pressure * $CONVERSIONTOPASCAL::T
  #############################
  @eval $fun_to(pressure::T) where T<:AbstractFloat = pressure * $CONVERSIONFROMPASCAL::T
  #############################
  @eval function $fun_from!(pressure::A)::Nothing where {T,A<:AbstractArray{<:T}}
  @simd for i in eachindex(pressure)
      @inbounds pressure[i] *= $CONVERSIONTOPASCAL
    end
  end
  ##############################
  @eval function $fun_to!(pressure::A)::Nothing where {T,A<:AbstractArray{<:T}}
    @simd for i in eachindex(pressure)
      @inbounds pressure[i] *= $CONVERSIONFROMPASCAL
    end
  end
end
