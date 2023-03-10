function G = calculate_sensitivity_indices(basis_index, coefficients, type)    

    N_p = size(basis_index, 2);

    norm = PCE_methods.get_multivariate_basis_norm(basis_index);

    c2_norm = sum(coefficients .* coefficients .* norm, 2);

    tot = sum(c2_norm);

    G = zeros(1, N_p);

    for k = 1:N_p

        if strcmp(type, 'total')
            filt = basis_index(:, k) > 0;
        else
            filt = sum(basis_index, 2) == basis_index(:, k);
        end

        c2_norm_p = c2_norm(filt);        
        G(k) = sum(c2_norm_p)/tot;
    end
end
