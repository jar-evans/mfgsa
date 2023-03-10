clear;
import PCE_methods.*; import general.*;

%% set time grid
T = 10;
dt = .01;
N_quad = int64(T/dt) + 1; %length(0:dt:T); % NOTE: should i change N_quad to
                                           % number of intervals and then
                                           % introduce new var N_points

%% generate realisations
input_ranges = [0.25 1.25 0.5];
input_means = [0.5 25/8 -1];

N_p = 3; grid_h = 0.2; % N = 1000;

[U, N] = general.generate_legendre_samples(grid_h, N_p);
% [U, N] = general.generate_legendre_samples(grid_h, N_p, N); 
I = general.generate_inputs(input_ranges, input_means, U);

%% generate N realisations of f, and then centre f
[fmk, t] = f_mechosc(I, T, dt);
mean_process = sum(fmk)/N;
fc = fmk - mean_process;  % should i be subtracting the mean process?
% fc = fmk - 0;  % OR the expectation of the mean process?

%% create projection matrix and calculate the PCE coeffs
N_ord = 4;
[coefficients, basis_index, proj_matrix] = PCE_methods.PCE(fc, U, N_ord);
% [basis_index, N_PC] = PCE_methods.generate_basis_index(N_ord, N_p);
% proj_matrix = PCE_methods.generate_projection_matrix(basis_index, U, N_p, N_PC, N_PC_quad);
% coefficients = PCE_methods.calculate_PCE_coefficients(proj_matrix, fc);


%% check approximations
y = PCE_methods.generate_approximation(proj_matrix, coefficients);
yp = PCE_methods.generate_approximation(proj_matrix, coefficients, sum(basis_index, 2) == 0);

figure(1);
[foo, bar] = f_mechosc(input_means, T, dt);
plot(t, y + mean_process); hold on;
plot(bar, foo); plot(t, yp + mean_process); plot(t, min(fmk));
plot(t, max(fmk)); hold off;

%% calculate sensitivity indices
G_tot_td = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, true, 'total', dt);
G_tot_pw = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, false, 'total', dt);

G_main_td = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, true, 'main', dt);
G_main_pw = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, false, 'main', dt);

figure(2);
plot(0:dt:T, G_tot_pw); hold on;
plot(0:dt:T, G_tot_td); hold off;

figure(3);
plot(0:dt:T, G_main_pw); hold on;
plot(0:dt:T, G_main_td); hold off;