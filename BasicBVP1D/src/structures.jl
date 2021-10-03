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

struct PeriodicBC
end
