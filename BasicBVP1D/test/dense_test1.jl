# solve the problem
a = 0;
b = 1;
n = 9;
f = 2*ones(n)
problem = FiniteDifferenceBVPProblem(a, b, n, f);
assemble_system!(problem);
u = solve_bvp(problem);

# computer error
err = @. problem.x * (1-problem.x) - u;
# check error 
norm(err, Inf)< 1e-14