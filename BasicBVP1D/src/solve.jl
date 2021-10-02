function solve_bvp(problem::FiniteDifferenceBVPProblem{TF, TBCL, TBCR}) where {TF, TBCL<:DirichletBC, TBCR<:DirichletBC}
    u = problem.A\problem.rhs;
    return u
end

function solve_bvp(problem::FiniteDifferenceBVPProblem{TF, TBCL, TBCR}) where {TF, TBCL<:NeumannBC, TBCR<:NeumannBC}
    u = lsmr(problem.A, problem.rhs);
    return u
end
