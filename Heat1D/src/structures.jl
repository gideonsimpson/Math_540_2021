struct FDHeatProblem{TBCL, TBCR,TTS}
    a # (a,b) define the domain
    b
    Δx # mesh spacing
    nx # number of points in the interior, excluding a and b
    x # mesh on which we solve the heat equation
    Δt # time step
    nΔt # number of time steps
    tmax # final time
    u0 # Initial condition
    A # the matrix (default to sparse)
    rhs # right hand side of our system
    left_bc::TBCL
    right_bc::TBCR
    time_step::TTS
end


struct DirichletBC
    u_bc
end

struct NeumannBC
    v_bc
end

struct PeriodicBC
end

struct BEuler
end