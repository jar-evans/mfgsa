clear; close all
addpath('../mfgsa')
addpath('../..\td_mfgsa\')

%% generate realisations
dt = 0.1; T = 10;
time_grid = 0:dt:T; N_t = length(time_grid);
N = 100; N_ord = 4;

input_ranges = [0.25 1.25 0.5]; input_means = [0.5 25/8 -1];
d = 3; 

U = general.generate_legendre_samples(d, N);
I = general.generate_inputs(input_ranges, input_means, U);

fmk = f_mechosc(I, time_grid);
mean_process = sum(fmk)/N;
fc = fmk - mean_process;

[coefficients, basis_index, proj_matrix] = PCE_methods.PCE(fc, U, N_ord);

% define high-fidelity and low-fidelity models
fcns{1} = @(Z, time_grid) f1(Z, time_grid);   % high-fidelity
% fcns{2} = @(Z) generate_realisations(mean_process, coefficients, basis_index, Z);
fcns{2} = @(Z, time_grid) f2(Z, time_grid);   % lowest-fidelity

w = [1; 0.1];%; 0.05; 0.001];       % assign model weights/costs
vec = [3 3];% 1 1];              % says each model is vectorized

budget = 2000;   % define computational budget

n_estimate = 100;
new_U = general.generate_legendre_samples(d, n_estimate);
Z = general.generate_inputs(input_ranges, input_means, new_U);
stats = estimate_statistics(fcns,Z,new_U,vec,time_grid);

%% COMPUTE MULTIFIDELITY GLOBAL SENSITIVITY ANALYSIS REPLICATES
n_reps = 100;    % number of replicates 
avg   = zeros(n_reps,N_t,2);    vr    = zeros(n_reps,N_t,2);
mc_sm = zeros(n_reps,N_t-1,d);    mc_st = zeros(n_reps,N_t-1,d);
mf_sm = zeros(n_reps,N_t-1,d);    mf_st = zeros(n_reps,N_t-1,d);

method = 'Owen';

for n = 1:n_reps

    n
    
    % call mfsobol.m with just the high-fidelity model to get Monte
    % Carlo estimate
    [sm,st,mu,sigsq] = mfsobol_vec(fcns(1),d,w(1),stats,budget,time_grid,vec(1),method);
    avg(n,:,1) = mu; 
    vr(n,:,1) = sigsq;
    mc_sm(n,:,:) = sm;
    mc_st(n,:,:) = st;
    
    % call mfsobol.m with full array of functions to get multifidelity
    % estimates
    [sm,st,mu,sigsq] = mfsobol_vec(fcns,d,w,stats,budget,time_grid,vec,method);
    avg(n,:,2) = mu; 
    vr(n,:,2) = sigsq;
    mf_sm(n,:,:) = sm;
    mf_st(n,:,:) = st;

end

%% PLOT ESTIMATOR SPREAD
warning('off','MATLAB:legend:IgnoringExtraEntries')
blue = [0       0.4470 0.7410];
red  = [0.8500  0.3250 0.0908];

% plot main effect sensitivity indices
figure(1); clf
h = boxplot([mc_sm(:,end,1), mf_sm(:,end,1), mc_sm(:,end,2), mf_sm(:,end,2), mc_sm(:,end,3), mf_sm(:,end,3)],...
    'Colors',[blue; red; blue; red; blue; red],'Whisker',10,...
    'labels',{'MC $s_m^1$','MF $s_m^1$','MC $s_m^2$','MF $s_m^2$','MC $s_m^3$','MF $s_m^3$'});
set(h,{'linew'},{2}); grid on
legend(flipud(findall(gca,'Tag','Box')), {'High-fidelity estimate','Multifidelity estimate'},...
    'Location','SouthWest','interpreter','latex'); legend boxoff
bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
title([method,' main sensitivity estimates for MO function'],'interpreter','latex')

if strcmp(method,'Owen') || strcmp(method,'Saltelli')
    % plot total effect sensitivity indices
    figure(2); clf
    h = boxplot([mc_st(:,end,1), mf_st(:,end,1), mc_st(:,end,2), mf_st(:,end,2), mc_st(:,end,3), mf_st(:,end,3)],...
        'Colors',[blue; red; blue; red; blue; red],'Whisker',10,...
        'labels',{'MC $s_t^1$','MF $s_t^1$','MC $s_t^2$','MF $s_t^2$','MC $s_t^3$','MF $s_t^3$'});
    set(h,{'linew'},{2}); grid on
    hLegend = legend(flipud(findall(gca,'Tag','Box')), {'High-fidelity estimate','Multifidelity estimate'},...
        'Location','NorthEast','interpreter','latex'); legend boxoff
    bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
    title([method,' total sensitivity estimates for MO function'],'interpreter','latex')
end

% plot variance estimates
% figure(3); clf
% histogram(vr(:,end,1),12,'facecolor',blue); hold on
% histogram(vr(:,end,2),12,'facecolor',red,'facealpha',1)
% legend({'High-fidelity','Multifidelity'},'interpreter','latex'); legend boxoff
% title('Variance estimates for Ishigami function','interpreter','latex')