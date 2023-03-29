function plot_variance(u, e, fc, v_PLOTS)

    if v_PLOTS
        var = sum(fc.*fc)/size(fc, 1);
        
        e2 = e.*e;
        var_KL_i = e2 .* transpose(u);
        var_KL = sum(var_KL_i, 2);
        
        figure(4);
        loglog(1:size(fc, 2), var); hold on;
        loglog(1:size(fc, 2), var_KL);
        xlim([1, size(fc, 2)]); hold off;
    end
end

