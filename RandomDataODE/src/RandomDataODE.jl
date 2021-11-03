module RandomDataODE
using DifferentialEquations

include("ode.jl")
export solve_ode
include("ensemble.jl")
export run_obs_ensemble, run_path_ensemble
end # module
