module BasicBVP1D

using SparseArrays
using LinearAlgebra
using FFTW

include("structures.jl")
export FiniteDifferenceBVPProblem, DirichletBC, NeumannBC, PeriodicBC, SpectralBVPProblem

include("constructors_fd.jl")
include("constructors_sparse.jl")
include("constructors_spectral.jl")
export FiniteDifferenceBVPProblem, SparseFiniteDifferenceBVPProblem, SpectralBVPProblem

include("assembly_fd.jl")
include("assembly_spectral.jl")
export assemble_system!

include("solve_fd.jl")
include("solve_spectral.jl")
export solve_bvp

end # module
