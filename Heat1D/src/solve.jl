function solve_heat(problem::FDHeatProblem{TBCL, TBCR, TS}) where {TBCL<:PeriodicBC, TBCR<:PeriodicBC, TS<:BEuler}
    u = similar(problem.u0);
    @. u = problem.u0;
    assemble_system!(problem);
    t = 0.;

    u_trajectory = [deepcopy(u)];
    t_trajectory = [t];

    for j in 1:problem.nΔt
        @. problem.rhs = u;
        u .= problem.A\problem.rhs;
        t+=problem.Δt;
        push!(u_trajectory, deepcopy(u));
        push!(t_trajectory, t);
    end

    return t_trajectory, u_trajectory
end