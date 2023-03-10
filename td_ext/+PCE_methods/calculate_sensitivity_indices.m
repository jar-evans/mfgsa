function G = calculate_sensitivity_indices(basis_index, coefficients, time_dependent, type, time_grid)

    if time_dependent
        G = generate_time_dep_indices(basis_index, coefficients, type, time_grid);
    else
        G = generate_pointwise_indices(basis_index, coefficients, type, time_grid);
    end

end

function G = generate_pointwise_indices(basis_index, coefficients, type, time_grid)

    N_p = size(basis_index, 2);
    N_quad = size(coefficients, 2);

    c2_tot = sum_over_bases(basis_index, coefficients, time_grid);

    tot = sum(c2_tot);

    G = zeros(N_p, N_quad);

    for k = 1:N_p

        if strcmp(type, 'total')

            filt = basis_index(:, k) > 0;
        else
            filt = sum(basis_index, 2) == basis_index(:, k);
        end

        c2_pk = c2_tot(filt, :);
        pk = sum(c2_pk);
        G(k, 1:N_quad) = pk./tot;
    end
end

function G = generate_time_dep_indices(basis_index, coefficients, type, time_grid)

    N_p = size(basis_index, 2);
    N_quad = size(coefficients, 2);

    c2_tot = sum_over_bases(basis_index, coefficients, time_grid);

    tot = sum(c2_tot);
    tot_cum = cumsum(tot);

    G = zeros(N_p, N_quad);

    for k = 1:N_p

        if strcmp(type, 'total')
            filt = basis_index(:, k) > 0;
        else
            filt = sum(basis_index, 2) == basis_index(:, k);
        end

        c2_pk = c2_tot(filt, :);
        pk = sum(c2_pk);
        pk_cum = cumsum(pk);
        G(k, 1:N_quad) = pk_cum./tot_cum;
    end
end
function c2_tot = sum_over_bases(basis_index, coefficients, time_grid)

    norm = PCE_methods.get_multivariate_basis_norm(basis_index);

    %w = dt; % uniform weighting

    w = max(time_grid)/length(time_grid);

    c2 = coefficients .* coefficients;
    c2 = w .* c2;
    c2_tot = c2 .* norm;

end


%             main_filt_1 = sum(basis_index, 2) == basis_index(:, 1);
%             main_filt_2 = sum(basis_index, 2) == basis_index(:, 2);
%             main_filt_3 = sum(basis_index, 2) == basis_index(:, 3);             
%             c2_p1_main = c2_tot(main_filt_1, :);
%             c2_p2_main = c2_tot(main_filt_2, :);
%             c2_p3_main = c2_tot(main_filt_3, :);
%             p1_main = sum(c2_p1_main);
%             p2_main = sum(c2_p2_main);
%             p3_main = sum(c2_p3_main);
%             p1_main_cum = cumsum(p1_main);
%             p2_main_cum = cumsum(p2_main);
%             p3_main_cum = cumsum(p3_main);
%             G1_main = p1_main_cum./tot_cum;
%             G2_main = p2_main_cum./tot_cum;
%             G3_main = p3_main_cum./tot_cum;
           

%             tot_filt_1 = basis_index(:, 1) > 0;
%             tot_filt_2 = basis_index(:, 2) > 0;
%             tot_filt_3 = basis_index(:, 3) > 0;
%             norm = PCE_methods.get_multivariate_basis_norm(basis_index);
%             w = dt; % uniform weighting
%             c2 = c .* c;
%             c2 = w .* c2;
%             c2_tot = c2 .* norm;
%             c2_p1 = c2_tot(tot_filt_1, :);
%             c2_p2 = c2_tot(tot_filt_2, :);
%             c2_p3 = c2_tot(tot_filt_3, :);
%             p1 = sum(c2_p1);
%             p2 = sum(c2_p2);
%             p3 = sum(c2_p3);
%             tot = sum(c2_tot);
%             p1_cum = cumsum(p1);
%             p2_cum = cumsum(p2);
%             p3_cum = cumsum(p3);
%             tot_cum = cumsum(tot);
%             G1 = p1_cum./tot_cum;
%             G2 = p2_cum./tot_cum;
%             G3 = p3_cum./tot_cum;

