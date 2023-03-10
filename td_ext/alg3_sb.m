clear;
import PCE_methods.*; import KL_methods.*; import general.*;

%% set time grid
T = 10;
dt = .1;
N_quad = int64(T/dt) + 1; %length(0:dt:T); % NOTE: should i change N_quad to
                                           % number of intervals and then
                                           % introduce new var N_points

%% turn plots on/off
u_PLOTS = false; 
e_PLOTS = false;
v_PLOTS = true;
approx_PLOTS = false;

%% generate realisations
input_ranges = [0.25 1.25 0.5];
input_means = [0.5 25/8 -1];

N_p = 3; grid_h = 0.2; % N = 1000;

[U, N] = general.generate_legendre_samples(grid_h, N_p);
% [U, N] = general.generate_legendre_samples(grid_h, N_p, N); 
I = general.generate_inputs(input_ranges, input_means, U);

%% generate N realisations of f, and then centre f
[fmk, t] = f(I, T, dt);
mean_process = sum(fmk)/N;
fc = fmk - mean_process;  % should i be subtracting the mean process?
% fc = fmk - 0;  % OR the expectation of the mean process?



%% form the covariance matrix (discretised covariance function) and compute dicretised KL modes
tol = 1e-1;
[u, e, N_KL] = KL_methods.eigen_stuff(N, N_quad, fc, tol);

KL_methods.plot_variance(u, e, fc, v_PLOTS);
KL_methods.plot_eigpairs(u, e, u_PLOTS, e_PLOTS);

fi = KL_methods.generate_KL_modes(fc, e);

approx = KL_methods.generate_approximation(fi, e, approx_PLOTS, mean_process);

N_ord = 4;
[coefficients, basis_index, ~] = KL_methods.PCE_KL_modes(fi, U, N_ord);

G_main = KL_methods.calculate_sensitivity_indices(basis_index, coefficients, 'main');
G_tot = KL_methods.calculate_sensitivity_indices(basis_index, coefficients, 'total');




















