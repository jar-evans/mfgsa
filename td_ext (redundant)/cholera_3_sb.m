clear;

T = 250;
dt = .05;%e-2;
% N_t = T/dt;

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

%% turn plots on/off
u_PLOTS = true; 
e_PLOTS = true;
v_PLOTS = true;
approx_PLOTS = true;

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

N_t = length(t);



%% form the covariance matrix (discretised covariance function) and compute dicretised KL modes
tol = 1e-6;
[u, e, N_KL] = KL_methods.eigen_stuff(N, N_t, fc, tol);

if v_PLOTS
    var = sum(fc.*fc)/size(fc, 1);

    e2 = e.*e;
    var_KL_i = e2 .* transpose(u);
    var_KL = sum(var_KL_i, 2);

    size(var_KL')
    size(var)

    figure(4);
    loglog(linspace(0, T, length(var)), var); hold on;
    loglog(linspace(0, T, length(var_KL)), var_KL');
    xlim([5e-2, 2e2]); hold off;
end

if (u_PLOTS)
    u = u;
    figure(1);
    plot(cumsum(u)/sum(u));

    figure(2);
    norm_u = u/max(u);
    semilogy(norm_u);
end

if (e_PLOTS)
    figure(3);
    semilogx(e);
end

fi = KL_methods.generate_KL_modes(fc, e);

% approx = KL_methods.generate_approximation(fi, e, approx_PLOTS, mean_process);

N_ord = 3;
[coefficients, basis_index, ~] = KL_methods.PCE_KL_modes(fi, U, N_ord);

G_main = KL_methods.calculate_sensitivity_indices(basis_index, coefficients, 'main');
G_tot = KL_methods.calculate_sensitivity_indices(basis_index, coefficients, 'total');




















