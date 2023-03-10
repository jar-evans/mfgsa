function G = gen_alg2(function_samples, time_grid, legendre_samples, PC_ord, time_dependent, type)
%% alg2
% This function uses 
%
%% Inputs
% *function_samples* - float^(N x N_k): output of model
% *time_grid* - 0:dt:T
% *legendre_samples* - float^(N x N_k): [-1, 1]
% *PC_ord* - int: order of PCE
% *time_dependent* - bool: true OR false i.e. generalised OR pointwise
% *type* - str: 'main' OR 'total' sensitivity indices
% 
%% Outputs
% *G* - float: sensitivity indices
% 
    [N, N_t] = size(function_samples);
    mean_process = sum(function_samples)/N;
    f_centred = function_samples - mean_process;  % should i be subtracting the mean process?
    % fc = function_samples - 0;  % OR the expectation of the mean process?

    %% create projection matrix and calculate the PCE coeffs
    [coefficients, basis_index, ~] = PCE_methods.PCE(f_centred, legendre_samples, PC_ord);

    %% calculate sensitivity indices
    G = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, time_dependent, type, time_grid);

end
