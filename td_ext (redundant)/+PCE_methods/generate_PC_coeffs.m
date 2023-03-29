function c = generate_PC_coeffs(fi, leg_samples, N_ord, N_p, basis_index, N_PC, N_PC_quad, N_KL, v)

    c = zeros(N_KL, N_PC);
    
    for i = 1:N_KL
        for n = 1:N_PC
            for j = 1:N_PC_quad

                c(i, n) = v*fi(j, i);

                for k = 1:N_p

                    res = legendre(basis_index(n, k), leg_samples(j, k), 'norm');

                    if (numel(res) > 1)
                        c(i, n) = c(i, n) * res(1);
                    else
                        c(i, n) = c(i, n) * res;
                    end

                end

            end
        end
    end

end