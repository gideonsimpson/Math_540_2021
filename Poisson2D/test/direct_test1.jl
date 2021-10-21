a = 0;
b = 1;
c = 0;
d = 1;
n = 9;

# set up exact solution and f
u_exact = (x,y) ->  x*(1-x)*y*(1-y)
f = (x,y)->2*y*(1-y)+2*x*(1-x)

# set up problem
problem = FDPoisson2DProblem(a,b,c,d,n,f);
# assemble system
assemble_system!(problem);
# solve
U = direct_solve_poisson2d(problem);
u = reshape(U,n,n);
# compute error
err = @. u - [u_exact(x_,y_) for x_ in problem.x, y_ in problem.y];
norm(err,Inf)<1e-14;
