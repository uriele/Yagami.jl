


#=

function raytracing!(results::RES,rays::RAY,refractive::M,knots_θ::V,knots_h::V,earthmodel::) where {T<:AbstractFloat, RES,RAY,  M<:AbstractMatrix{T}, V<:AbstractVector{T}}

    # Initialize results
    @inbounds for i in eachindex(results)
        results[i] = ResultRay{T}()
    end

    # Loop through rays
    @inbounds for i in eachindex(rays)
        ray = rays[i]
        point_x = __getpointx(ray)
        point_y = __getpointy(ray)
        direction_x = __getdirectionx(ray)
        direction_y = __getdirectiony(ray)

        # Perform intersection calculations
        # This is a placeholder for the actual intersection logic
        # You would replace this with your fast marching method or other intersection logic
        intersection_length = 0.0  # Placeholder for length calculation
        intersection_altitude = 0.0  # Placeholder for altitude calculation
        intersection_azimuth = 0.0  # Placeholder for azimuth calculation

        # Store results
        results[i] = ResultRay{T}(direction_x, direction_y, point_x, point_y,
                                   indices_i=0, indices_j=0, length=intersection_length,
                                   altitude=intersection_altitude, azimuth=intersection_azimuth)
    end

    return results  # Return the populated results array or vector of intersec

end
=#
