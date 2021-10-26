function serial_pi_estimation(n_samples)
    n_inside = 0;

    x = zeros(2);

    for _ in 1:n_samples
        @. x = rand();
        if(norm(x,2)<1)
            n_inside+=1;
        end
    end

    pi_estimate = 4 * n_inside/n_samples;

    return pi_estimate

end

function pi_estimation(n_samples)
    n_inside = 0;

    x_inside = zeros(Bool, n_samples);

    x = zeros(2);
    Threads.@threads for j in 1:n_samples
        @. x = rand();
        if(norm(x,2)<1)
            x_inside[j] = 1;
        end
    end

    n_inside = sum(x_inside);

    pi_estimate = 4 * n_inside/n_samples;

    return pi_estimate

end

function atomic_pi_estimation(n_samples)

    n_inside = Atomic{Int}(0);

    x = zeros(2);

    Threads.@threads for _ in 1:n_samples
        @. x = rand();
        if(norm(x,2)<1)
            atomic_add!(n_inside,1);
        end
    end

    pi_estimate = 4 * n_inside[]/n_samples;

    return pi_estimate

end