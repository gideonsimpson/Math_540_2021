# solve the problem
a = 0;
b = 1;
n = 7;
f = x-> (π)^2 * sin(π * x);
u_exact = x-> sin(π * x);

problem = SpectralBVPProblem(a, b, n, f, DirichletBC(0), DirichletBC(0));
assemble_system!(problem);
u = solve_bvp(problem);

# computer error
err = @. u - u_exact(problem.x);
# check error 
norm(err, Inf)< 1e-14