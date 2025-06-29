
@inline function snellslaw(normalx::T,normaly::T, directionx::T,directiony::T, n_incident::T=T(1.0),n_transmitted::T=T(1.0)) where {T<:AbstractFloat}
  if n_incident==n_transmitted # no bending
     return directionx, directiony
  end
  n01 = n_incident/n_transmitted
  n01² = n01*n01
  cosθ_incident = -clamp(directionx*normalx+directiony*normaly, -1, 1)
  sinθ²_transmitted = n01²*(1-cosθ_incident*cosθ_incident)
  if sinθ²_transmitted ≤ 1
    tmp = (n01*cosθ_incident - sqrt(1-sinθ²_transmitted))
    return (n01*directionx + tmp*normalx), (n01*directiony + tmp*normaly)
  else
    return (directionx - 2*cosθ_incident*normalx), (directiony - 2*cosθ_incident*normaly)
  end
end

function snellslaw!(::Type{T}=Float64;normal::V,direction::V,n_incident::T=T(1.0),n_transmitted::T=T(1.0)) where {T<:AbstractFloat,V<:AbstractVector{T}}
  return __bend_ray(n_incident,n_transmitted,normal,direction)
end

@inline function __bend_ray(n_incident::T,n_transmitted::T,normal::V,direction::V) where {T<:AbstractFloat,V<:AbstractVector{T}}
  local n01=n_incident/n_transmitted
  local n01²=n01*n01
  local cosθ_incident=-clamp(dot(direction,normal),-1,1)
  local sinθ²_transmitted =n01²*(1-cosθ_incident*cosθ_incident)

  if sinθ²_transmitted ≤ 1
    tmp=(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))
    direction[1]= n01*direction[1]+tmp*normal[1]
    direction[2]= n01*direction[2]+tmp*normal[2]
  else
    direction[1]= direction[1]-2*cosθ_incident*normal[1]
    direction[2]= direction[2]-2*cosθ_incident*normal[2]
  end
  return true
end
