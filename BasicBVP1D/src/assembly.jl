"""
`assemble_system!` - assemble the linear system for solving
"""
function assemble_system!(problem::FiniteDifferenceBVPProblem{TF, TBCL, TBCR}) where {TF<:AbstractArray, TBCL<:DirichletBC, TBCR<:DirichletBC}

    # build first row
    problem.A[1,1] = 2/problem.Δx^2;
    problem.A[1,2] = -1/problem.Δx^2;
    # build rows 2-n-1
    for i in 2:problem.n-1
        problem.A[i,i-1] = -1/problem.Δx^2;
        problem.A[i,i] = 2/problem.Δx^2;
        problem.A[i,i+1] = -1/problem.Δx^2;
    end
    # build row n
    problem.A[problem.n,problem.n-1] =  -1/problem.Δx^2;
    problem.A[problem.n,problem.n] =  2/problem.Δx^2;

    # build rhs
    @. problem.rhs = problem.f;
    problem.rhs[1] += problem.left_bc.u_bc/problem.Δx^2
    problem.rhs[problem.n] += problem.right_bc.u_bc/problem.Δx^2

    problem
end

function assemble_system!(problem::FiniteDifferenceBVPProblem{TF, TBCL, TBCR}) where {TF<:Function, TBCL<:DirichletBC, TBCR<:DirichletBC}

    # build first row
    problem.A[1,1] = 2/problem.Δx^2;
    problem.A[1,2] = -1/problem.Δx^2;
    # build rows 2-n-1
    for i in 2:problem.n-1
        problem.A[i,i-1] = -1/problem.Δx^2;
        problem.A[i,i] = 2/problem.Δx^2;
        problem.A[i,i+1] = -1/problem.Δx^2;
    end
    # build row n
    problem.A[problem.n,problem.n-1] =  -1/problem.Δx^2;
    problem.A[problem.n,problem.n] =  2/problem.Δx^2;

    # build rhs
    @. problem.rhs = problem.f(problem.x);
    problem.rhs[1] += problem.left_bc.u_bc/problem.Δx^2
    problem.rhs[problem.n] += problem.right_bc.u_bc/problem.Δx^2

    problem
end
