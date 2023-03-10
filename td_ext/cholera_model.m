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
[t, y] = f_cholera(I, dt, T, ic, tols, N, 2);


figure(1);
loglog(0:dt:T, y);
ylim([5e-2 1e4]);
xlim([5e-2 2e2]);

% figure(1);
% loglog(t, y(:, 1:3));
% ylim([5e-2 1e4]);
% xlim([5e-2 2e2]);
% 
% figure(2);
% loglog(t, y(:, 4:5));
% xlim([5e-2 2e2]);

% %% set time grid % time discr may not always be unique
% N_t = length(time_grid);    % #points in time discretisation
% N_h = N_t - 1;               % #intervals in time discretisation
% 
% G_alg2 = gen_alg2(fmk, time_grid, U, 4, true, 'main');
% 
% G_alg3 = gen_alg3(fmk, time_grid, U, 1e-5, 4, 'main');

G_alg3 = gen_alg3(y, t, U, 1e-3, 3, 'total')
% figure(2);
% semilogx(0:dt:T, G_alg2);

