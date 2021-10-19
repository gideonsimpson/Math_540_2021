function direct_solve_poisson2d(problem::FDPoisson2DProblem) 
    u = problem.A\problem.rhs;
    return u
end
