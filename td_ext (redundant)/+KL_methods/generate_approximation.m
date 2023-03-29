function approx = generate_approximation(modes, e, approx_PLOTS, mean_process)

    [N_quad, N_KL] = size(e); 
    N = size(modes, 1);
    approx = zeros(N, N_quad);

    %modes = transpose(modes);

    for i = 1:N_KL
        approx = approx + modes(:, i)*transpose(e(:, i));
    end

    approx = sum(approx);

    mp = false;
    if exist('mean_process', 'var')
        approx = approx + mean_process;
        mp = true;
    end

    if approx_PLOTS
        plot(1:N_quad, approx, 'go'); hold on;
        if mp
            plot(1:N_quad, mean_process); 
        end
        hold off;
    end

end

