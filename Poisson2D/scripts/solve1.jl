using Poisson2D
using JLD2
using Printf

# set problem parameters
a = 0;
b = 1;
c = 0;
d = 1;
n = 9;

kx = 5;
ky = 4;
f = (x,y)->((kx*π/(b-a))^2 + (ky * π/(d-c))^2) * sin(kx*π*(x-a)/(b-a)) * sin(ky*π*(y-c)/(d-c));

# set up and solve problem
problem = FDPoisson2DProblem(a,b,c,d,n,f);
assemble_system!(problem);
U = direct_solve_poisson2d(problem);

# name of file
fname = @sprintf("soln1_n%d.jld2", n);

# save to disk
jldsave(fname; problem, U);