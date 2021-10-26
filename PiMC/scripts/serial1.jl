using PiMC
using BenchmarkTools
using Random
using Printf

Random.seed!(100);

n_samples = 10^8;

Random.seed!(100);
pi_est = serial_pi_estimation(n_samples);
@printf("π estimate with %.1e samples = %g\n", n_samples,pi_est);

@btime serial_pi_estimation(n_samples);