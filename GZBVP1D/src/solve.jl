function solve_gz(u0, problem::FDGZProblem)

    function F!(F, u)
        assemble_system!(u, problem);
        @. F = problem.F;
        F
    end

    soln = nlsolve(F!, u0);
    return soln
end