function proj_matrix = generate_projection_matrix(basis_index, U)
    
    [N, N_p] = size(U);
    N_PC = size(basis_index, 1);
    
    v = N_p/N;
    proj_matrix = zeros(N_PC, N);
    
    for n = 1:N_PC
        for j = 1:N
            proj_matrix(n, j) = v;
            for k = 1:N_p
                norm = 2/((2*basis_index(n, k)) + 1);
                res = legendre(basis_index(n, k), U(j, k))/norm;
    
                if (numel(res) > 1)
                    proj_matrix(n, j) = proj_matrix(n, j) * res(1);
                else
                    proj_matrix(n, j) = proj_matrix(n, j) * res;
                end
    
            end
    
        end
    end

end