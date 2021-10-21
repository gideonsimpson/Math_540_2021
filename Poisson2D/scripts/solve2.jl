using Poisson2D
using LinearAlgebra
using BenchmarkTools

@show n_threads = Threads.nthreads();
# without the next line, it will use the max number of threads available
BLAS.set_num_threads(n_threads);
# set problem parameters
a = 0;
b = 1;
c = 0;
d = 1;
n = 999;

kx = 5;
ky = 4;
f = (x,y)->((kx*π/(b-a))^2 + (ky * π/(d-c))^2) * sin(kx*π*(x-a)/(b-a)) * sin(ky*π*(y-c)/(d-c));

# set up and solve problem
problem = FDPoisson2DProblem(a,b,c,d,n,f);
assemble_system!(problem);
@btime U = direct_solve_poisson2d(problem);
