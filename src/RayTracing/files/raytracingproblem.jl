# Available data formats for ray tracing problems
# Add new formats here as needed.
abstract type AvailableDataFormat end
struct NCFormat <: AvailableDataFormat end
struct GeofitFormat <: AvailableDataFormat end

const FORMATFILES=(:ncformat,:geofitformat)

# Add new data formats here as needed.
# Examples:
# struct NewDataFormat <: AvailableDataFormat end
#
# Create the file newdataformat.jl in the src/RayTracing/files directory
# add the name as a symbol to FORMATFILES (e.g. :newdataformat)
##################################################################################
# List of available data formats
# add to the list as new formats are implemented
# at the end of the file, this value will be used to include the necessary files
# the naming convention in use is <nameofthedataformat>.jl


"""
  $SIGNATURES

Data structure for the ray tracing problem.
This struct holds all the necessary data for performing ray tracing calculations, including the atmosphere, refractive index, and satellite scan information.
"""
struct RayTracingProblem{T<:AbstractFloat,MT<:MeanType,AM<:AirModel,
  EA<:EarthApproximation,ATM<:AtmosphereSetting,
  V<:AbstractVector{T},
  M<:AbstractMatrix{T},L<:AbstractLogger
}
  filename::String
  meantype::MT
  model::AM
  earthmodel::EA
  ###### Data
  atmosphere::ATM
  refractive::M
  pointsx::V
  pointsy::V
  directionsx::V
  directionsy::V
  ##############
  tangent_h:: V
  tangent_θ:: V
  nscans::Int
  nlos::Int
  __logger::L
end

"""
   $SIGNATURES
Create a `RayTracingProblem` from an available data format. This function is a factory method that dispatch
 a `RayTracingProblem` based on the specified data format.

 # Arguments:
 - `fileorfolder::String`: The path to the file or folder containing the ray tracing data.
  - `dataformat::ADF`: The data format of the file or folder, which can be `NCFormat` or `MipasFormat`. (default is `NCFormat`).
 - `meantype::MT=GeometricMean()`: The mean type to be used for the refractive index calculation (default is `GeometricMean`).
 - `model::AM=Ciddor()`: The air model to be used for the refractive index calculation (default is `Ciddor`).
  - `earthmodel::EA=Fukushima()`: The earth approximation model to be used (default is `Fukushima`).
  - `logger=nothing`: An optional logger to log information during the problem creation (default is `NullLogger`).

  # Returns:
  - `RayTracingProblem`: An instance of `RayTracingProblem` containing the data read from the file or folder and the calculated refractive index and atmosphere.
"""
RayTracingProblem(fileorfolder::String,dataformat::ADF=NCFormat();
meantype::MT=GeometricMean(), model::AM=Ciddor(),
earthmodel::EA=Fukushima(),
logger=nothing) where {MT<:MeanType,AM<:AirModel,
EA<:EarthApproximation,ADF<:AvailableDataFormat} = throw(ArgumentError("Unknown data format: $dataformat. Use NCFormat or MipasFormat."))



for format in FORMATFILES
  include("$(format).jl")
end
##### Interfaces for different data formats #####
RayTracingProblem(fileorfolder::String,::NCFormat;kwargs...)    = NCRayTracingProblem(fileorfolder;kwargs...)
RayTracingProblem(fileorfolder::String,::MipasFormat;kwargs...) = GeofitRayTracingProblem(fileorfolder;kwargs...)
#################################################################################################################
