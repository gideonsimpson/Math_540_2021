module BasicBVP1D

using SparseArrays
using LinearAlgebra

include("structures.jl")
export FiniteDifferenceBVPProblem, DirichletBC, NeumannBC

include("constructors.jl")
export FiniteDifferenceBVPProblem, SparseFiniteDifferenceBVPProblem

include("assembly.jl")
export assemble_system!

include("solve.jl")
export solve_bvp

end # module
