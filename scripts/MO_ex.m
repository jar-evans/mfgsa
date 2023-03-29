%% SETUP
clear; close all
addpath('../mfgsa')
addpath('../td_ext\')

%% generate realisations
dt = 0.1; T = 10; time_grid = 0:dt:T; N_t = length(time_grid);
input_ranges = [0.25 1.25 0.5]; input_means = [0.5 25/8 -1];

budget = 1e3; d = 3;


%% COMPUTE MULTIFIDELITY GLOBAL SENSITIVITY ANALYSIS REPLICATES
n_reps = 50;    % number of replicates 

mc_sm = zeros(n_reps,N_t,d);    mc_st = zeros(n_reps,N_t,d);

method = 'Owen';

fcns{1} = @(Z) f1(Z, time_grid);   % high-fidelity


p = budget/(d+2);
for n = 1:n_reps

    Z_A = generate_inputs(p);
    U_A = norm_samples(Z_A);

    Z_B = generate_inputs(p);
    U_B = norm_samples(Z_B);

    yA = fcns{1}(Z_A);
    yB = fcns{1}(Z_B);

    yC = zeros(p,N_t,d);
    for i = 1:d
        Z_Ci = Z_B;
        Z_Ci(:,i) = Z_A(:,i);
        yC(:,:,i) = fcns{1}(Z_Ci);
    end

    [sm,st] = estimate_sobol(method,yA,yB,yC,true,true);

    mc_sm(n,:,:) = sm;
    mc_st(n,:,:) = st;

    n


end


mc_sm = reshape(mean(mc_sm, 1), N_t, d);
mc_sm = mc_sm(end, :)

mc_st = reshape(mean(mc_st, 1), N_t, d);
mc_st = mc_st(end, :)