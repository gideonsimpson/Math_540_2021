function serial_pi_estimation(n_samples)
    n_inside = 0;

    for j in 1:n_samples
        x = rand(2);
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

    Threads.@threads for j in 1:n_samples
        x = rand(2);
        if(norm(x,2)<1)
            x_inside[j] = 1;
        end
    end

    n_inside = sum(x_inside);

    pi_estimate = 4 * n_inside/n_samples;

    return pi_estimate

end