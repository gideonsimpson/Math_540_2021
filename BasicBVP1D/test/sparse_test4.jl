f = x-> 0;
u_exact = x-> x;

a = 0;
b = 1;
n = 9;

problem = SparseFiniteDifferenceBVPProblem(0, 1, n, f, NeumannBC(1), NeumannBC(1));
assemble_system!(problem);
u = solve_bvp(problem);
norm(@. u - problem.x - u[1] )<1e-12