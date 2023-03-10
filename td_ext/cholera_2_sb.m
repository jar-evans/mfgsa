clear;

T = 250;
dt = .05;%e-2;
N_quad = T/dt;

N_pop = 10000;
S0 = N_pop - 1; I0 = 1; R0 = 0;
B_H0 = 0; B_L0 = 0;
tol_abs = 1e-6;
tol_rel = 1e-6;

% uncertain params
bL = 1.5;
bH = 7.5;
kL = 1e6;
kH = kL/700;
b = 1/1560;
X = 168/5;
Z = 70;
d = 7/30;
g = 7/5;

ic = [N_pop, S0, I0, R0, B_H0, B_L0];
tols = [tol_abs, tol_rel];


%% generate realisations

input_means = [bL, bH, kL, kH, b, X, Z, d, g];
input_ranges = 0.2.*input_means;

N_p = length(input_means); grid_h = 0; N = 50;
% 
% [U, N] = general.generate_legendre_samples(grid_h, N_p);
[U, N] = general.generate_legendre_samples(grid_h, N_p, N); 
I = general.generate_inputs(input_ranges, input_means, U);

%% generate N realisations of f, and then centre f
[t, fmk] = f_cholera(I, dt, T, ic, tols, N, 2);
mean_process = sum(fmk)/N;
fc = fmk - mean_process;  % should i be subtracting the mean process?


%% create projection matrix and calculate the PCE coeffs
N_ord = 3;
[coefficients, basis_index, proj_matrix] = PCE_methods.PCE(fc, U, N_ord);
% [basis_index, N_PC] = PCE_methods.generate_basis_index(N_ord, N_p);
% proj_matrix = PCE_methods.generate_projection_matrix(basis_index, U, N_p, N_PC, N_PC_quad);
% coefficients = PCE_methods.calculate_PCE_coefficients(proj_matrix, fc);


%% check approximations
y = PCE_methods.generate_approximation(proj_matrix, coefficients);
yp = PCE_methods.generate_approximation(proj_matrix, coefficients, sum(basis_index, 2) == 0);

figure(1);
[foo, bar] = f_cholera(mean(I), dt, T, ic, tols, 1, 2);
plot(t, y + mean_process); hold on;
plot(foo, bar); plot(t, yp + mean_process); plot(t, min(fmk));
plot(t, max(fmk)); hold off;

%% calculate sensitivity indices
G_tot_td = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, true, 'total', dt);
G_tot_pw = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, false, 'total', dt);

G_main_td = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, true, 'main', dt);
G_main_pw = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, false, 'main', dt);

% figure(2);
% semilogx(0:dt:T, G_tot_pw);

figure(4);
semilogx(0:dt:T, G_tot_td);

% figure(3);
% semilogx(0:dt:T, G_main_pw); hold on;
% semilogx(0:dt:T, G_main_td); hold off;