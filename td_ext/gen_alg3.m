function G = gen_alg3(function_samples, time_grid, legendre_samples, KL_tol, PC_ord, type)
%% alg3
% This function uses 
%
%% Inputs
% *function_samples* - float^(N x N_k): output of model
% *time_grid* - 0:dt:T
% *legendre_samples* - float^(N x N_k): [-1, 1]
% *KL_tol* - float: determines truncation level of KL expansion
% *PC_ord* - int: order of PCE
% *type* - str: 'main' OR 'total' sensitivity indices
% 
%% Outputs
% *G* - float: sensitivity indices
% 
    [N, N_t] = size(function_samples);
    mean_process = sum(function_samples)/N;
    f_centred = function_samples - mean_process;  % should i be subtracting the mean process?
    % fc = function_samples - 0;  % OR the expectation of the mean process?

    %% form the covariance matrix (discretised covariance function) and compute dicretised KL modes
    [~, e, ~] = KL_methods.eigen_stuff(N, N_t, f_centred, KL_tol);
    fi = KL_methods.generate_KL_modes(f_centred, e);
    [coefficients, basis_index, ~] = KL_methods.PCE_KL_modes(fi, legendre_samples, PC_ord);
    
    G = KL_methods.calculate_sensitivity_indices(basis_index, coefficients, type);
end