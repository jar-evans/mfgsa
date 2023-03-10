function plot_eigpairs(u, e, u_PLOTS, e_PLOTS)
    if (u_PLOTS)
        figure(1);
        plot(cumsum(u)/sum(u));
        
        figure(2);
        norm_u = u/max(u);
        plot(log(norm_u));
    end
    
    if (e_PLOTS)
        figure(3);
        plot(e); 
    end

end

