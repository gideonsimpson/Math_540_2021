using Poisson2D
using BenchmarkTools

@show Threads.nthreads();
# set problem parameters
a = 0;
b = 1;
c = 0;
d = 1;
n = 99;

kx = 5;
ky = 4;
f = (x,y)->((kx*π/(b-a))^2 + (ky * π/(d-c))^2) * sin(kx*π*(x-a)/(b-a)) * sin(ky*π*(y-c)/(d-c));

# set up and solve problem
problem = FDPoisson2DProblem(a,b,c,d,n,f);
assemble_system!(problem);
@btime U = direct_solve_poisson2d(problem);
