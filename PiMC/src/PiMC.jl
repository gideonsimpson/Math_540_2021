module PiMC

using LinearAlgebra
using Base.Threads

include("estimators.jl")
export serial_pi_estimation, pi_estimation, atomic_pi_estimation

end # module
