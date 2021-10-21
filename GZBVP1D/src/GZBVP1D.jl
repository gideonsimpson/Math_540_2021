module GZBVP1D

using LinearAlgebra
using NLsolve

include("structures.jl")
export FDGZProblem, DirichletBC, NeumannBC, PeriodicBC

include("constructors.jl")
export FDGZProblem

include("assembly.jl")
export assemble_system!

include("solve.jl")
export solve_gz

end # module
