
# right hand side function in our system
function f!(du, u, p, t)
    du[1] = u[1] * (1.0-u[2]);
    du[2] = u[2] * (u[1]-1.0);
    du
end

# set up ODE
function solve_ode(u0, Tmax)
    # define the ODE problem
    prob = ODEProblem(f!, u0, (0., Tmax))
    # solve the problem
    sol = solve(prob)

    return sol
end