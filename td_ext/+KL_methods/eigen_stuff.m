function [u, e, N_KL] = eigen_stuff(N, N_quad, fc, tol)

    Klm = zeros(N_quad);
    
    for k = 1:N               %optimise
        for l = 1:N_quad
            for m = 1:N_quad
                Klm(l, m) = Klm(l, m) + fc(k,l)*fc(k,m);      
            end
        end
    end
    
    Klm = Klm/(N-1); % R^N_quad x N_quad
    
    %% solve the eigenpair problem (W^1/2 K W^1/2 ui = lambda ui)
    %W = diag(w);
    
    [e, U] = eig(Klm);  % [vectors, values] (for w = 1)
    [u, vec_order] = sort(diag(U), 'descend');
    ind = 1:length(u);
    
    %N_KL = N_quad; % non truncated form
    norm_u = u/max(u);
    ind = ind(norm_u > tol);
    vec_order = vec_order(norm_u > tol);
    
    u = u(ind);
    e = e(:,vec_order);
    N_KL = length(u); % truncated form


end