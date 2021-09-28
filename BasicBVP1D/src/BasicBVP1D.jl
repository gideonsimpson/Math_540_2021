module BasicBVP1D

using SparseArrays
struct FiniteDifferenceBVPProblem{TF}
    a # (a,b) define the domain
    b
    Δx # mesh spacing
    n # number of points in the interior, excluding a and b
    x # interior mesh
    f::TF # value of f at the interior nodes
    A # the matrix
    rhs # right hand side of our system
end

"""
`FiniteDifferenceBVPProblem` - Initalize the data structure
"""
function FiniteDifferenceBVPProblem(a, b, n, f)
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2)[2:end-1];

    # allocate the matrix
    A = zeros(n, n);
    rhs = zeros(n);
    
    return FiniteDifferenceBVPProblem(a, b, Δx, n, x, f, A, rhs)
end


"""
`SparseFiniteDifferenceBVPProblem` - Initalize the data structure
"""
function SparseFiniteDifferenceBVPProblem(a, b, n, f)
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2)[2:end-1];

    # allocate the matrix
    A = spzeros(n, n);
    rhs = zeros(n);
    
    return FiniteDifferenceBVPProblem(a, b, Δx, n, x, f, A, rhs)
end


"""
`assemble_system!` - assemble the linear system for solving
"""
function assemble_system!(problem::FiniteDifferenceBVPProblem{TF}) where {TF<:AbstractArray}

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

    problem
end

function assemble_system!(problem::FiniteDifferenceBVPProblem{TF}) where {TF<:Function}

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

    problem
end


function solve_bvp(problem)
    u = problem.A\problem.rhs;
    return u
end

export FiniteDifferenceBVPProblem, assemble_system!, solve_bvp, SparseFiniteDifferenceBVPProblem
end # module
