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
