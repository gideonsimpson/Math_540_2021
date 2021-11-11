using Random
using Distributions
using BenchmarkTools
using RandomDataODE

n_samples = 10^4;
Tmax = 6.;

# generate initial condition data
ϵ = 0.2;
xdist = Uniform(0.5-ϵ, 0.5 + ϵ);
ydist = Uniform(2.0-ϵ, 2.0 + ϵ);
Random.seed!(100)
u0data = [[rand(xdist), rand(ydist)] for _ in 1:n_samples];

f = sol->sol(Tmax)[1]; # observable function


@btime samples = run_obs_ensemble(u0data, Tmax, f);