clear; close all
addpath('../mfgsa')
addpath('../..\td_mfgsa\')

%% generate realisations
%% set grid
dt = 0.001; T = 250; time_grid = 0:dt:T;
N_t = length(time_grid); N_quad = N_t -1;

% uncertain params
bH = 7.5;
kL = 1e6;
X = 168/5;
Z = 70;
g = 7/5;

%% set sample settings
input_mean = [bH, kL, X, Z, g]; d = length(input_mean);
input_range = 0.2.*input_mean;

N_p = length(input_mean);

% U = general.generate_legendre_samples(d, N);
% I = general.generate_inputs(input_ranges, input_means, U);

% define high-fidelity and low-fidelity models
fcns{1} = @(I, time_grid) f_cholera(I, 0:dt:T, 2, 4);   % high-fidelity
fcns{2} = @(I, time_grid) f_cholera(I, 0:dt:T, 2, 1); % lowest-fidelity   % lowest-fidelity

w = [1; 0.25];%; 0.05; 0.001];       % assign model weights/costs
vec = [3 3];% 1 1];              % says each model is vectorized

budget = 100;   % define computational budget

n_estimate = 10;
new_U = general.generate_legendre_samples(d, n_estimate);
Z = general.generate_inputs(input_range, input_mean, new_U);
stats = estimate_statistics(fcns,Z,new_U,vec,time_grid);

%% COMPUTE MULTIFIDELITY GLOBAL SENSITIVITY ANALYSIS REPLICATES
n_reps = 50;    % number of replicates 
avg   = zeros(n_reps,N_t,2);    vr    = zeros(n_reps,N_t,2);
mc_sm = zeros(n_reps,N_t-1,d);    mc_st = zeros(n_reps,N_t-1,d);
mf_sm = zeros(n_reps,N_t-1,d);    mf_st = zeros(n_reps,N_t-1,d);

method = 'Owen';

for n = 1:n_reps

    n
    
    % call mfsobol.m with just the high-fidelity model to get Monte
    % Carlo estimate
    [sm,st,mu,sigsq] = mfsobol_vec(fcns(1),d,w(1),stats,budget,time_grid,vec(1),method,input_mean,input_range);
    avg(n,:,1) = mu; 
    vr(n,:,1) = sigsq;
    mc_sm(n,:,:) = sm;
    mc_st(n,:,:) = st;
    
    % call mfsobol.m with full array of functions to get multifidelity
    % estimates
    [sm,st,mu,sigsq] = mfsobol_vec(fcns,d,w,stats,budget,time_grid,vec,method,input_mean,input_range);
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
title([method,' main sensitivity estimates for SIR function'],'interpreter','latex')

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
    title([method,' total sensitivity estimates for SIR function'],'interpreter','latex')
end

% plot variance estimates
% figure(3); clf
% histogram(vr(:,end,1),12,'facecolor',blue); hold on
% histogram(vr(:,end,2),12,'facecolor',red,'facealpha',1)
% legend({'High-fidelity','Multifidelity'},'interpreter','latex'); legend boxoff
% title('Variance estimates for SIR function','interpreter','latex')