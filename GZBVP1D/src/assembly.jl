function assemble_system!(u, problem::FDGZProblem{TBC, TBC}) where {TBC<:DirichletBC}


    problem.F[1] =problem.ϵ * (u[2] - 2*u[1] + problem.left_bc.u_bc)/problem.Δx^2 +  u[1] * (1-u[1]^2);

    for i in 2:problem.n-1
        problem.F[i]= problem.ϵ * (u[i+1] - 2*u[i] +u[i-1])/problem.Δx^2 +u[i] * (1-u[i]^2);
    end
    problem.F[end] = problem.ϵ * (problem.right_bc.u_bc - 2*u[end] + u[end-1])/problem.Δx^2 +  u[end] * (1-u[end]^2);

    problem
end
