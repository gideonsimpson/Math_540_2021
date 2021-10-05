"""
`SpectralBVPProblem` - Initalize the data structure
"""
function SpectralBVPProblem(a, b, n, f, bcl::TBCL, bcr::TBCR) where {TBCL<:PeriodicBC, TBCR<:PeriodicBC}
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2)[1:end-1];
    N = n+1;

    # allocate the matrix
    λ = zeros(n+1);
    rhs = zeros(n+1);

    return SpectralBVPProblem(a, b, Δx, n, N, x, f, λ, rhs, bcl, bcr)
end

function SpectralBVPProblem(a, b, n, f, bcl::TBCL, bcr::TBCR) where {TBCL<:DirichletBC, TBCR<:DirichletBC}
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2)[2:end-1];
    N = 2*(n+1);

    # allocate the matrix
    λ = zeros(N);
    rhs = zeros(N);

    return SpectralBVPProblem(a, b, Δx, n, N, x, f, λ, rhs, bcl, bcr)
end

function SpectralBVPProblem(a, b, n, f, bcl::TBCL, bcr::TBCR) where {TBCL<:NeumannBC, TBCR<:NeumannBC}
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2);
    N = 2*(n+1);

    # allocate the matrix
    λ = zeros(N);
    rhs = zeros(N);

    return SpectralBVPProblem(a, b, Δx, n, N, x, f, λ, rhs, bcl, bcr)
end