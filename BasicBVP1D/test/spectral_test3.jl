# solve the problem
a = 0;
b = 1;
n = 7;
f = x-> (π)^2 * cos(π * x);
u_exact = x-> cos(π * x);

problem = SpectralBVPProblem(a, b, n, f, NeumannBC(0), NeumannBC(0));
assemble_system!(problem);
u = solve_bvp(problem);

# computer error
err = @. u - u_exact(problem.x);
# check error 
norm(err, Inf)< 1e-14