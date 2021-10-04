function solve_bvp(problem::SpectralBVPProblem{TF, TBCL, TBCR}) where {TF, TBCL<:PeriodicBC, TBCR<:PeriodicBC}
    rhshat = fft(problem.rhs);
    uhat = similar(rhshat);
    @. uhat[2:end] = rhshat[2:end]/problem.λ[2:end];
    u = ifft(uhat);
    return u
end

function solve_bvp(problem::SpectralBVPProblem{TF, TBCL, TBCR}) where {TF, TBCL<:DirichletBC, TBCR<:DirichletBC}
    rhshat = fft(problem.rhs);
    uhat = similar(rhshat);
    @. uhat[2:end] = rhshat[2:end]/problem.λ[2:end];
    u_ = ifft(uhat)
    u = u_[problem.N÷2 + 2:end];
    return u
end

function solve_bvp(problem::SpectralBVPProblem{TF, TBCL, TBCR}) where {TF, TBCL<:NeumannBC, TBCR<:NeumannBC}
    rhshat = fft(problem.rhs);
    uhat = similar(rhshat);
    @. uhat[2:end] = rhshat[2:end]/problem.λ[2:end];
    u_ = ifft(uhat);
    u = [u_[problem.N÷2 + 1:end]; u_[1]];
    return u
end

