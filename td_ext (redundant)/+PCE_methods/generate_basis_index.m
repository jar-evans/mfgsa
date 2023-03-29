function [basis_index, N_PC] = generate_basis_index(N_ord, N_p)
    N_PC = factorial(N_ord + N_p)/(factorial(N_ord)*factorial(N_p));
    basis_index = unique(nchoosek(repmat(0:N_ord, 1, N_p), N_p), 'rows');
    basis_index = basis_index(sum(basis_index, 2) <= N_ord, :);
end






