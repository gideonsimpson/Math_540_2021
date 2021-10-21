struct FDGZProblem{TBCL, TBCR}
    a # (a,b) define the domain
    b
    Δx # mesh spacing
    n # number of points in the interior, excluding a and b
    x # interior mesh
    ϵ # scalar parameter
    F # holds the residual values
    left_bc::TBCL
    right_bc::TBCR
end


struct DirichletBC
    u_bc
end

struct NeumannBC
    v_bc
end

struct PeriodicBC
end
