function [coefficients, basis_index, proj_matrix] = PCE(evaluations, U, N_ord, basis_index, proj_matrix)

    N_p = size(U, 2);

    if ~exist('basis_index', 'var')
        basis_index = PCE_methods.generate_basis_index(N_ord, N_p);
    end

    if ~exist('proj_matrix', 'var')
        proj_matrix = PCE_methods.generate_projection_matrix(basis_index, U);
    end

    coefficients = PCE_methods.calculate_PCE_coefficients(proj_matrix, evaluations);




%     [basis_index, N_PC] = PCE_methods.generate_basis_index(N_ord, N_p);
%     
%     N_KL = size(fi, 1);
%     
%     % A_range = 0.25;
%     % B_range = 1.25;
%     % L_range = 0.5;
%     
%     v = N_p/N_PC_quad;%1;%A_range*B_range*L_range/N_PC_quad;
%     c = PCE_methods.generate_PC_coeffs(fi, leg_samples, N_ord, N_p, basis_index, N_PC, N_PC_quad, N_KL, v);

end


% N_ord = 4; N_p = 3; N_PC_quad = N;
% [basis_index, N_PC] = PCE_methods.generate_basis_index(N_ord, N_p);
% proj_matrix = PCE_methods.generate_projection_matrix(basis_index, U, N_p, N_PC, N_PC_quad);
% coefficients = PCE_methods.calculate_PCE_coefficients(proj_matrix, fc);