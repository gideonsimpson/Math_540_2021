
"""
`SparseFiniteDifferenceBVPProblem` - Initalize the data structure
"""
function SparseFiniteDifferenceBVPProblem(a, b, n, f)  
    return SparseFiniteDifferenceBVPProblem(a, b, n, f, DirichletBC(0), DirichletBC(0))
end

function SparseFiniteDifferenceBVPProblem(a, b, n, f, bcl::TBCL, bcr::TBCR) where {TBCL<:DirichletBC, TBCR<:DirichletBC}
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2)[2:end-1];

    # allocate the matrix
    A = spzeros(n, n);
    rhs = zeros(n);
    
    return FiniteDifferenceBVPProblem(a, b, Δx, n, x, f, A, rhs,
        bcl, bcr)
end

function SparseFiniteDifferenceBVPProblem(a, b, n, f, bcl::TBCL, bcr::TBCR) where {TBCL<:NeumannBC, TBCR<:NeumannBC}
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2);

    # allocate the matrix
    A = spzeros(n+2, n+2);
    rhs = zeros(n+2);
    
    return FiniteDifferenceBVPProblem(a, b, Δx, n, x, f, A, rhs,
        bcl, bcr)
end

function SparseFiniteDifferenceBVPProblem(a, b, n, f, bcl::TBCL, bcr::TBCR) where {TBCL<:PeriodicBC, TBCR<:PeriodicBC}
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2)[1:end-1];

    # allocate the matrix
    A = spzeros(n+1, n+1);
    rhs = zeros(n+1);
    
    return FiniteDifferenceBVPProblem(a, b, Δx, n, x, f, A, rhs, 
        bcl, bcr)
end
