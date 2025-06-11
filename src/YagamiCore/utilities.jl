for datuminfo in DATUMINFO
  wgsdatum = Symbol("WGS84", datuminfo)
  refdatum = Symbol("REF", datuminfo)
  funcdatum = datuminfo

 @eval $funcdatum(::Type{<:T}=Float64) where T<:AbstractFloat = convert(T, $refdatum[])
 @eval $funcdatum(::Type{<:Real})  = convert(T, $refdatum[])
end


# setdatum! function to set the datum parameters
@eval function setdatum!(datum::Datum=WGS84Latest)::Nothing where Datum
    if datum ∉ subtypes(Datum)
      @warn("Invalid datum type: $datum. Using current parameterization.")
      return nothing
    end
     _EARTH = ellipsoid(datum)

    _MAJORAXIS = majoraxis(_EARTH) |> er-> uconvert(km,er) |> ustrip
    _MINORAXIS = minoraxis(_EARTH) |> er-> uconvert(km,er) |> ustrip
    _ECCENTRICITY² = eccentricity²(_EARTH) |> ustrip
    _ECCENTRICITY = sqrt(_ECCENTRICITY²)   |> ustrip
    _NORMALIZEMINORAXIS = _MINORAXIS / _MAJORAXIS
    _COMPLECCENTRICITY² = 1 - _ECCENTRICITY²
    $(SETDATUM...)
    return nothing
end

# setmajoraxis! function to set the major axis
@eval function setmajoraxis!(majoraxis::T=REFMAJORAXIS[])::Nothing where T<:AbstractFloat
      _MAJORAXIS = majoraxis
      _MINORAXIS = REFMINORAXIS[] # use the current minor axis
      if _MAJORAXIS < _MINORAXIS
        @warn("Major axis must be greater than or equal to minor axis. But got: $majoraxis < $(REFMINORAXIS[]). Keeping current values ($REFMAJORAXIS, $REFMINORAXIS).")
        return nothing
      end
      _ECCENTRICITY² = 1 - (_MINORAXIS / _MAJORAXIS)^2
      _ECCENTRICITY = sqrt(_ECCENTRICITY²)
      _NORMALIZEMINORAXIS = _MINORAXIS / _MAJORAXIS
      _COMPLECCENTRICITY² = 1 - _ECCENTRICITY²
      $(SETDATUM...)
    return nothing
end

# setminoraxis! function to set the minor axis
@eval function setminoraxis!(minoraxis::T=REFMINORAXIS[])::Nothing where T<:AbstractFloat
      _MINORAXIS = minoraxis
      _MAJORAXIS = REFMAJORAXIS[] # use the current major axis
      if _MINORAXIS > _MAJORAXIS
        @warn("Minor axis must be less than or equal to major axis. But got: $minoraxis > $(REFMAJORAXIS[]). Keeping current values ($REFMAJORAXIS, $REFMINORAXIS).")
        return nothing
      end
      _ECCENTRICITY² = 1 - (_MINORAXIS / _MAJORAXIS)^2
      _ECCENTRICITY = sqrt(_ECCENTRICITY²)
      _NORMALIZEMINORAXIS = _MINORAXIS / _MAJORAXIS
      _COMPLECCENTRICITY² = 1 - _ECCENTRICITY²
      $(SETDATUM...)
    return nothing
end

# seteccentricity²! function to set the eccentricity squared
@eval function seteccentricity²!(eccentricity²::T=REFECCENTRICITY²[])::Nothing where T<:AbstractFloat
      _ECCENTRICITY² = eccentricity²
      if _ECCENTRICITY² < 0 || _ECCENTRICITY² > 1
        @warn("Eccentricity squared must be in the range [0, 1]. But got: $eccentricity². Keeping current values ($(REFECCENTRICITY²[])).")
        return nothing
      end
      _ECCENTRICITY = sqrt(_ECCENTRICITY²)
      _MAJORAXIS = REFMAJORAXIS[] # use the current major axis
      _MINORAXIS = _MAJORAXIS * sqrt(1 - _ECCENTRICITY²)
      _NORMALIZEMINORAXIS = _MINORAXIS / _MAJORAXIS
      _COMPLECCENTRICITY² = 1 - _ECCENTRICITY²
      $(SETDATUM...)
    return nothing
end
