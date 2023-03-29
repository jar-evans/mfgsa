function y = f2(I, t)

    a = I(:, 1);
    b = I(:, 2);
    l = I(:, 3);

    y = l.*exp(-a*t).*(cos(1.1*b*t) + 0.8*(a/b)*sin(b*t));

% 
%     %like if this is a surrogate
% %     it should be 
% 
% i = 0
% if i < var
% 
% T = dt = 
% time_grid = 
% 
% ...
%     ...
%     ...
%     #
% 
% coeff, basis etc = PCE
% 
% i += 1
% 
% if i < var
% 
%     gen realisations (c, b, etc)
% 
% %% set time grid
% T = 10;
% dt = .1;
% N_quad = int64(T/dt) + 1
% 
% %% generate realisations
% input_ranges = [0.25 1.25 0.5];
% input_means = [0.5 25/8 -1];
% 
% N_p = 3; grid_h = 0;  N = 100;
% 
% [U, N] = general.generate_legendre_samples(grid_h, N_p, N); 
% I = general.generate_inputs(input_ranges, input_means, U);
% 
% %% generate N realisations of f, and then centre f
% [fmk, t] = f_mechosc(I, T, dt);
% mean_process = sum(fmk)/N;
% fc = fmk - mean_process;
% %% create projection matrix and calculate the PCE coeffs
% N_ord = 4;
% [coefficients, basis_index, proj_matrix] = PCE_methods.PCE(fc, U, N_ord);
% [new_U, ~] = general.generate_legendre_samples(0, 3, 100);
% generate_realisations(mean_process, coefficients, basis_index, U);
% 
% 
% %% calculate sensitivity indices
% G_tot_td = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, true, 'total', dt);
% G_tot_pw = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, false, 'total', dt);
% 
% G_main_td = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, true, 'main', dt);
% G_main_pw = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, false, 'main', dt);
% 
