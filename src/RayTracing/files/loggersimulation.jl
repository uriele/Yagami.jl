using UnPack:@unpack
using Logging,LoggingExtras


struct DebugLogger{L<:AbstractLogger} <: AbstractLogger
  logger::L
end

# Forwarded methods for DebugLogger
Logging.min_enabled_level(logger::DebugLogger) =
    min_enabled_level(logger.logger)

Logging.shouldlog(logger::DebugLogger, level, mod, group, id) =
    shouldlog(logger.logger, level, mod, group, id)

Logging.handle_message(logger::DebugLogger, level, message, mod, group, id, file, line) =
    handle_message(logger.logger, level, message, mod, group, id, file, line)

Logging.with_logger(f::Function, logger::DebugLogger) =
    with_logger(f, logger.logger)

Logging.global_logger(logger::DebugLogger) =
    global_logger(logger.logger)

function DebugLogger(filename::String)
  DebugLogger(infologger(filename))
end


@inline write_line_logger(logger::L,zb::ZB,ray::Int,iter::Int) where {L<:AbstractLogger, ZB<:Zbrent} = nothing

@inline function write_line_logger(logger::DL,zb::ZB,ray::Int,iter::Int) where {DL<:DebugLogger, ZB<:Zbrent}
  with_logger(logger) do
    iter==1 && @info "Ray: $ray"
    i= zb.f.i
    j= zb.f.j
    pointx = zb.f.pointx
    pointy = zb.f.pointy
    directionx = zb.f.directionx
    directiony = zb.f.directiony
    altitude,azimuth=__innerfunc(zb,pointx,pointy)
    fmin= zb.fx
    length_t   = zb.x
    a   = zb.a
    b   = zb.b
    niter= zb.__iter

    spacing=textshort(" ";padding=4)
    rowstr = [
            spacing,
            numshort(iter,        padding=10),
            numshort(i,           padding=8),
            numshort(j,           padding=8),
            numshortf(pointx,      padding=8),
            numshortf(pointy,      padding=8),
            numshortf(directionx,  padding=8),
            numshortf(directiony,  padding=8),
            numshortf(altitude,    padding=8),
            numshortf(azimuth,     padding=8),
            numshortf(length_t,    padding=8),
            numshortf(a, padding=10),
            numshortf(b, padding=10),
            numshort(fmin, padding=10),
            numshort(niter, padding=10),
    ]
    @info join(rowstr, "")
  end
end



@inline function write_tracing_log(logger::L, results::RR,refractive::M) where {L<:AbstractLogger,
  T<:AbstractFloat,RR<:AbstractMatrix{<:AbstractResult{T}}, M<:AbstractMatrix{T}
}

  Ni,Mi=size(refractive)
  with_logger(logger) do
    @info " Completed Ray Tracing Simulation"
    @info ""
    @info "============================================================"
    @info "Ray Tracing Results"
    @info "============================================================"
    spacing=textshort(" ";padding=4)
    header = [
      spacing,
      textshort("iter";          padding=10),
      textshort("i";             padding=8),
      textshort("j";             padding=8),
      textshort("x [km]";        padding=8),
      textshort("y [km]";        padding=8),
      textshort("dX";            padding=8),
      textshort("dY";            padding=8),
      textshort("h [km]";        padding=8),
      textshort("θ [°]";         padding=8),
      textshort("length [km]";   padding=8),
      textshort("(n-1)";         padding=10),
      textshort("islevel";       padding=10),
      textshort("isdescending";  padding=12),
    ]


    @info join(header, "")
    for ray in axes(results, 2) # iterate over the rays
      @info "Ray: $(ray)/$(size(results,2))"
      for iter in axes(results, 1) # iterate over the iterations

        @unpack i,j,pointx,pointy,directionx,directiony,altitude,azimuth,length_t,islevel,descending = results[iter,ray]

        ni1= FREESPACE
        if (1<=i<=Ni) && (1<=j<=Mi)
          ni1 = refractive[i,j]
        end
        islvl = islevel ? "level" : "radius"
        isdesc = descending ? "descending" : "ascending"
        rowstr = [
          spacing,
          numshort(iter,        padding=10),
          numshort(i,           padding=8),
          numshort(j,           padding=8),
          numshortf(pointx,      padding=8),
          numshortf(pointy,      padding=8),
          numshortf(directionx,  padding=8),
          numshortf(directiony,  padding=8),
          numshortf(altitude,    padding=8),
          numshortf(azimuth,     padding=8),
          numshortf(length_t,    padding=8),
          numshort((ni1-1), padding=10),
          textshort(islvl,    padding=10),
          textshort(isdesc,     padding=12),
        ]
        @info join(rowstr, "")
        # stop when you reach the end of the rays
        if iter>1 && length_t == 0
          break # break the loop if the length is zero
        end
      end
    end
    @info "============================================================"
  end
end
@inline write_tracing_log(::NullLogger, ::RR,::M) where {T<:AbstractFloat,RR<:AbstractMatrix{<:AbstractResult{T}},M<:AbstractMatrix{T}} = nothing

@inline write_tracing_log(::DL, ::RR,::M) where {T<:AbstractFloat,
  RR<:AbstractMatrix{<:AbstractResult{T}},
  M<:AbstractMatrix{T}, DL<:DebugLogger} = nothing




function atmosphere_info(atmosphere)
      @info "===================================================================="
      @info "Temperature information:"
      @info "===================================================================="
      @info "$(textshort("θ [°]")) $(textshort("h [km]")) $(textshort("T [K]"))  $(textshort("P [Pa]")) $(textshort("humidity [%]")) $(textshort("CO2 [ppm]")) $(textshort("λ [μm]"))"
      for j in eachindex(atmosphere.temperature.knots_h)
        h= atmosphere.temperature.knots_h[j]
        for i in eachindex(atmosphere.temperature.knots_θ[1:end-1])
          θ= atmosphere.temperature.knots_θ[i]
          T= atmosphere.temperature(θ,h)
          P= atmosphere.pressure(θ,h)
          H= atmosphere.humidity(θ,h)*100
          CO2= atmosphere.co2ppm(θ,h)
          λ= atmosphere.wavelength(θ,h)
          @info "$(numshort(θ)) $(numshort(h)) $(numshort(T)) $(numshort(P)) $(numshort(H)) $(numshort(CO2)) $(numshort(λ))"
        end
      end
end

function refractive_info(refractive)
  @info "===================================================================="
  @info "Refractive index information:"
  @info "===================================================================="
  @info "$(textshort("θ_i")) $(textshort("h_j")) $(textshort("(n-1) [1e-4]"))"
  for j in axes(refractive,2)
    for i in axes(refractive,1)
      n= refractive[i,j]
      @info "$(numshort(i)) $(numshort(j)) $(numshort((n-1)*1e4))"
    end
  end
end
