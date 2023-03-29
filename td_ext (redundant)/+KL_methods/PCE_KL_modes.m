function [coefficients, basis_index, proj_matrix] = PCE_KL_modes(fi, U, N_ord)

    N_KL = size(fi, 1);

    [ci, basis_index, proj_matrix] =  PCE_methods.PCE(transpose(fi(1, :)), U, N_ord);
    N_PC = size(basis_index, 1);

    coefficients = zeros(N_PC, N_KL);
    coefficients(:, 1) = ci;

    for i = 2:N_KL
        [ci, ~, ~] = PCE_methods.PCE(transpose(fi(i, :)), U, N_ord, basis_index, proj_matrix);
        coefficients(:, i) = ci;
    end

end
