function solve_bvp(problem)
    u = problem.A\problem.rhs;
    return u
end
