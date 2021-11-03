

# compute an ensemble of ϕ(u) scalar observable values
function run_obs_ensemble(u0data, Tmax, ϕ)

    n = length(u0data);
    samples = zeros(n);

    Threads.@threads for i in 1:n
        soln = solve_ode(u0data[i], Tmax);
        samples[i] = ϕ(soln)
    end

    return samples
end


# compute an ensemble of the paths at the indicated values of 
function run_path_ensemble(u0data, Tmax, t_vals)

    n = length(u0data);
    u1_samples = zeros(n, length(t_vals));
    u2_samples = zeros(n, length(t_vals));

    Threads.@threads for i in 1:n
        soln = solve_ode(u0data[i], Tmax);
        soln_at_t = soln(t_vals);
        u1_samples[i,:] .= [u[1] for u in soln_at_t];
        u2_samples[i,:] .= [u[2] for u in soln_at_t];        
    end

    return u1_samples,u2_samples
end

