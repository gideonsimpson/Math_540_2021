
function assemble_system!(problem::FDHeatProblem{TBCL, TBCR, TS}) where {TBCL<:PeriodicBC, TBCR<:PeriodicBC, TS<:BEuler}

    # build first row
    problem.A[1,1] = 1 + 2*problem.Δt/problem.Δx^2;
    problem.A[1,2] = -problem.Δt/problem.Δx^2;
    problem.A[1,end] = -problem.Δt/problem.Δx^2;
    # build rows 2-n
    for i in 2:problem.nx
        problem.A[i,i-1] = -problem.Δt/problem.Δx^2;
        problem.A[i,i] = 1 + 2*problem.Δt/problem.Δx^2;
        problem.A[i,i+1] = -problem.Δt/problem.Δx^2;
    end
    # build row nx+1
    problem.A[end,1] =  -problem.Δt/problem.Δx^2;
    problem.A[end,end-1] =  -problem.Δt/problem.Δx^2;
    problem.A[end,end] =  1 + 2*problem.Δt/problem.Δx^2;

    problem
end
