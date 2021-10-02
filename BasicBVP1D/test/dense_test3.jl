f = x->0;
u_exact = x->x;

a = 0;
b = 1;
n = 9;

problem = FiniteDifferenceBVPProblem(0, 1, n, f, DirichletBC(0), DirichletBC(1));
assemble_system!(problem);
u = solve_bvp(problem);
err = @. problem.x - u;
# check error 
norm(err, Inf)< 1e-14