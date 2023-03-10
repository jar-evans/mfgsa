function norm = get_multivariate_basis_norm(basis_index)
    
    [N_PC, N_p] = size(basis_index);
    norm = ones(N_PC, 1);

    for i = 1:N_PC
        for k = 1:N_p
            n = basis_index(i, k);
            norm(i) = norm(i) * (2./(2.*n + 1));
        end
    end
end