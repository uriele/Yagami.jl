module CurtisGodson
  using DocStringExtensions
  using LinearAlgebra,StructArrays,FastGaussQuadrature
  include("quadraturepoints.jl")
  export QuadraturePoints, get_quadrature_points, linintegral

end
