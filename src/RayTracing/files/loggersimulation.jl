using UnPack:@unpack
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
write_tracing_log(::NullLogger, ::RR,::M) where {T<:AbstractFloat,RR<:AbstractMatrix{<:AbstractResult{T}},M<:AbstractMatrix{T}} = nothing
