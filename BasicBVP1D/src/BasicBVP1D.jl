module BasicBVP1D

struct FiniteDifferenceBVPProblem
    a # (a,b) define the domain
    b
    Δx # mesh spacing
    n # number of points in the interior, excluding a and b
    x # interior mesh
    f # value of f at the interior nodes
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

export FiniteDifferenceBVPProblem
end # module
