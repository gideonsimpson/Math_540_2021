module Poisson2D

using LinearAlgebra
using SparseArrays

include("structures.jl")
export FDPoisson2DProblem, DirichletBC, NeumannBC, PeriodicBC
include("constructors_sparse.jl")
export FDPoisson2DProblem
include("assembly_fd.jl")
export assemble_system!
include("solve_fd.jl")
export direct_solve_poisson2d

end # module
