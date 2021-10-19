struct FDPoisson2DProblem{TF, TBCT, TBCB, TBCL, TBCR}
    a # (a,b)×(c,d) define the domain
    b
    c
    d
    Δx # mesh spacing in x
    nx # number of points in the interior, excluding a and b
    Δy # mesh spacing in y
    ny # number of points in the interior, excluding c and d
    x # interior mesh
    y # interior mesh
    f::TF # value of f at the interior nodes
    A # the matrix
    rhs # right hand side of our system
    top_bc::TBCT
    bottom_bc::TBCB
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
