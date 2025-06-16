using Logging,LoggingExtras
using Test
using Yagami.RayTracing
logfile="$(pwd())/RayTracing/_data/cairt.log"
testfile="$(pwd())/RayTracing/_data/cairt.nc"
@info "Testing Cairt File: $testfile"
problem=RayTracingProblem(testfile;logger=nothing)

@test getlogger(problem) isa NullLogger
problem=RayTracingProblem(testfile;logger=logfile)
@test getlogger(problem) isa FormatLogger

#@testset "Test Cairt Extrapolation"
  @test num_edges(problem)[1]==length(problem.temperature.knots_θ)-1
  @test num_edges(problem)[2]==length(problem.temperature.knots_h)
  for totest in (:temperature, :pressure, :humidity, :co2ppm,:wavelength)
      extfunc = Symbol("extrapolate", totest)
      nodefunc = Symbol("node_", totest)
      @testset "Testing helper $totest" begin
        @eval begin
          @test problem.$(totest)(1.0,2.0) == $extfunc(problem, 1.0, 2.0)
          for j in 1:num_edges(problem)[2]
            for i in 1:num_edges(problem)[1]
              local knot_hj = problem.temperature.knots_h[j]
              local knot_θi = problem.temperature.knots_θ[i]
              p1=problem.$(totest)(knot_θi,knot_hj)
              p2=$nodefunc(problem, i, j)[2]
              if isnan(p1)
                @test isnan(p2)
              else
                @test p1 == p2
              end
            end
          end
        end
    end
  end

# note refractive has an extra col for wrapping around without mod1
  @testset "Testing Wedge Refractive Index" begin
  @test num_wedges(problem) == size(problem.refractive[1:end-1,:])

  for j in 1:num_wedges(problem)[2]
    for i in 1:num_wedges(problem)[1]
      n1=wedge_refractive(problem,i,j)[2]
      n2=problem.refractive[i,j]
      if isnan(n1)
        @test isnan(n2)
      else
        @test n1 == n2
      end
    end
  end
end
