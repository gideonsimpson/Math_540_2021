# solve the problem
a = 0;
b = 1;
n = 7;
f = x-> (2*π)^2 * sin(2*π * x);
u_exact = x-> sin(2*π * x);

problem = SpectralBVPProblem(a, b, n, f, PeriodicBC(), PeriodicBC());
assemble_system!(problem);
u = solve_bvp(problem);
# computer error
err = @. u - u_exact(problem.x);
# check error 
norm(err, Inf)< 1e-14