%% SETUP
clear; close all
addpath('../mfgsa')
addpath('../../td_mfgsa\')

%% generate realisations
dt = 0.05; T = 250; time_grid = 0:dt:T; N_t = length(time_grid);
 % uncertain params
bH = 7.5;
kL = 1e6;
X = 168/5;
Z = 70;
g = 7/5;

%% generate realisations

input_mean = [bH, kL, X, Z, g];
input_range = 0.2.*input_mean;

N_p = length(input_mean);

budget = 1e5; d = 5;
 N = budget;

%% COMPUTE MULTIFIDELITY GLOBAL SENSITIVITY ANALYSIS REPLICATES
n_reps = 2;    % number of replicates 

mc_sm = zeros(n_reps,d);    mc_st = zeros(n_reps,d);

method = 'Owen';

fcns{1} = @(I) f_cholera(I, time_grid, 2, 45);   % high-fidelity

f = @(I) f_cholera(I, time_grid, 2, 45);


p = floor(budget/(d+2));
for n = 1:n_reps

    U_A = general.generate_legendre_samples(N_p, p); 
    Z_A = general.generate_model_inputs(input_range, input_mean, U_A);

    U_B = general.generate_legendre_samples(N_p, p); 
    Z_B = general.generate_model_inputs(input_range, input_mean, U_B);

    yA = fcns{1}(Z_A);
    yB = fcns{1}(Z_B);

    yC = zeros(p,N_t,d);
    for i = 1:d
        Z_Ci = Z_B;
        Z_Ci(:,i) = Z_A(:,i);
        yC(:,:,i) = fcns{1}(Z_Ci);
    end

    [sm,st] = estimate_sobol(method,yA,yB,yC,true,true);

    mc_sm(n,:) = sm(end,:)
    mc_st(n,:) = st(end,:)

%     base = TD_SURROGATE(f, time_grid, input_range, input_mean); N_ord = 2; N_KL = 15;
% 
%         N = budget;
% 
%         kl_surrogate = base.generate_KLPC_surrogate(N, N_KL, N_ord);
%         G_tot_kl = kl_surrogate.calculate_sensitivity_indices('total')
%         G_main_kl = kl_surrogate.calculate_sensitivity_indices('main');
% 
%         a = sum(kl_surrogate.eigenvalues)




end


mc_sm = mean(mc_sm, 1);

mc_st = mean(mc_st, 1);