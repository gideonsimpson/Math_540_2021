function assemble_system!(problem::SpectralBVPProblem{TF, TBCL, TBCR}) where {TF<:Function, TBCL<:PeriodicBC, TBCR<:PeriodicBC}

    # build eigenvalues
    @. problem.λ = ((2*π/(problem.b-problem.a)) * [0:problem.N÷2; -(problem.N÷2)+1:-1])^2;

    # build rhs
    @. problem.rhs = problem.f(problem.x);

    problem
end

function assemble_system!(problem::SpectralBVPProblem{TF, TBCL, TBCR}) where {TF<:AbstractArray, TBCL<:DirichletBC, TBCR<:DirichletBC}

    # build eigenvalues
    @. problem.λ = ((π/(problem.b-problem.a)) * [0:problem.N÷2; -(problem.N÷2)+1:-1])^2;
    
    # build rhs
    @. problem.rhs = problem.f;

    problem
end

function assemble_system!(problem::SpectralBVPProblem{TF, TBCL, TBCR}) where {TF<:Function, TBCL<:DirichletBC, TBCR<:DirichletBC}

    # build eigenvalues
    @. problem.λ = ((π/(problem.b-problem.a)) * [0:problem.N÷2; -(problem.N÷2)+1:-1])^2;

    # build rhs with odd extension
    @. problem.rhs[problem.N÷2+1:end] = problem.f([problem.a; problem.x]);
    problem.rhs[1:problem.N÷2] .= -reverse(problem.f.([problem.x; problem.b]));

    problem
end


function assemble_system!(problem::SpectralBVPProblem{TF, TBCL, TBCR}) where {TF<:AbstractArray, TBCL<:NeumannBC, TBCR<:NeumannBC}

    # build eigenvalues
    @. problem.λ = ((π/(problem.b-problem.a)) * [0:problem.N÷2; -(problem.N÷2)+1:-1])^2;
    
    # build rhs
    @. problem.rhs = problem.f;

    problem
end

function assemble_system!(problem::SpectralBVPProblem{TF, TBCL, TBCR}) where {TF<:Function, TBCL<:NeumannBC, TBCR<:NeumannBC}

    # build eigenvalues
    @. problem.λ = ((π/(problem.b-problem.a)) * [0:problem.N÷2; -(problem.N÷2)+1:-1])^2;

    # build rhs with even extension
    @. problem.rhs[problem.N÷2+1:end] = problem.f(problem.x[1:end-1]);
    problem.rhs[1:problem.N÷2] .= reverse(problem.f.(problem.x[2:end]));

    problem
end

