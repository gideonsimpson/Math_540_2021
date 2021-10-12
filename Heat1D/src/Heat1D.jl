module Heat1D

using LinearAlgebra
using SparseArrays

include("structures.jl")
export FDHeatProblem, PeriodicBC, BEuler
include("constructors.jl")
export FDHeatProblem
include("assembly.jl")
export assemble_system!
include("solve.jl")
export solve_heat

end # module
