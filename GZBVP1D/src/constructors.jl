"""
`FDGZProblem` - Initalize the data structure
"""
function FDGZProblem(ϵ, n)
    return FDGZProblem(0, 1, ϵ, n,  DirichletBC(-1), DirichletBC(1))
end

function FDGZProblem(a, b, ϵ, n, bcl::TBCL, bcr::TBCR) where {TBCL<:DirichletBC, TBCR<:DirichletBC}
    Δx = (b-a)/(n+1);
    x = LinRange(a, b, n+2)[2:end-1];
    F = zeros(n);
    
    return FDGZProblem(a, b, Δx, n, x, ϵ, F, bcl, bcr)
end
