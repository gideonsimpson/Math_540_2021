module BasicBVP1D

using SparseArrays
struct FiniteDifferenceBVPProblem{TF, TBCL, TBCR}
    a # (a,b) define the domain
    b
    Δx # mesh spacing
    n # number of points in the interior, excluding a and b
    x # interior mesh
    f::TF # value of f at the interior nodes
    A # the matrix
    rhs # right hand side of our system
    left_bc::TBCL
    right_bc::TBCR
end


struct DirichletBC
    u_bc
end
struct NeumannBC
    v_bc
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
    
    return FiniteDifferenceBVPProblem(a, b, Δx, n, x, f, A, rhs, 
        DirichletBC(0), DirichletBC(0))
end

function FiniteDifferenceBVPProblem(a, b, n, f, bcl::TBCL, bcr::TBCR) where {TBCL<:DirichletBC, TBCR<:DirichletBC}
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2)[2:end-1];

    # allocate the matrix
    A = zeros(n, n);
    rhs = zeros(n);
    
    return FiniteDifferenceBVPProblem(a, b, Δx, n, x, f, A, rhs, 
        bcl, bcr)
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
    
    return FiniteDifferenceBVPProblem(a, b, Δx, n, x, f, A, rhs,
        DirichletBC(0), DirichletBC(0))
end

"""
`SparseFiniteDifferenceBVPProblem` - Initalize the data structure
"""
function SparseFiniteDifferenceBVPProblem(a, b, n, f, bcl::TBCL, bcr::TBCR) where {TBCL<:DirichletBC, TBCR<:DirichletBC}
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2)[2:end-1];

    # allocate the matrix
    A = spzeros(n, n);
    rhs = zeros(n);
    
    return FiniteDifferenceBVPProblem(a, b, Δx, n, x, f, A, rhs,
        bcl, bcr)
end


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


function solve_bvp(problem)
    u = problem.A\problem.rhs;
    return u
end

export FiniteDifferenceBVPProblem, assemble_system!, solve_bvp, SparseFiniteDifferenceBVPProblem, DirichletBC, NeumannBC
end # module
