function w = generate_realisations(mean_process, coefficients, basis_index, inputs)

    [N_PC, N_t] = size(coefficients);
    [N, N_p] = size(inputs);
    w = zeros(N, N_t);

    for j = 1:N
        for n = 1:N_PC
           
            res = ones(1,N_t);

            for k = 1:N_p
                norm = 2/((2*basis_index(n, k)) + 1);
                temp = legendre(basis_index(n, k), inputs(j, k));
                res = res .* temp(1,:) ./ norm;
            end
            w(j,:) = sum([w(j,:); res .* coefficients(N_PC,:)]);
        end
    end

    w = w + mean_process;
    
end
