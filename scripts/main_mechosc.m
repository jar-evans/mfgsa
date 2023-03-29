clear; close all
addpath('../mfgsa')
addpath('../td_ext\')

%% generate realisations
dt = 0.1; T = 10;
time_grid = 0:dt:T; N_t = length(time_grid);

input_ranges = [0.25 1.25 0.5]; input_means = [0.5 25/8 -1];

d = 3; N = 1000; N_ord = 4;

w = [1, 0.005];%; 0.05; 0.001];       % assign model weights/costs
vec = [1, 3];% 1 1];              % says each model is vectorized

budget = 200;   % define computational budget

m = [budget-N,N];

%% COMPUTE MULTIFIDELITY GLOBAL SENSITIVITY ANALYSIS REPLICATES
n_reps = 100;    % number of replicates 
avg   = zeros(n_reps,N_t,2);    vr    = zeros(n_reps,N_t,2);
mc_sm = zeros(n_reps,N_t,d);    mc_st = zeros(n_reps,N_t,d);
mf_sm = zeros(n_reps,N_t,d);    mf_st = zeros(n_reps,N_t,d);

mf_smb = mf_sm;

method = 'Owen';

for n = 1:n_reps
% 
    [U, N] = general.generate_legendre_samples(0, d, N);
    I = general.generate_inputs(input_ranges, input_means, U);

    fmk = f_mechosc(I, time_grid);
    mean_process = sum(fmk)/N;
    fc = fmk - mean_process;

    [coefficients, basis_index, ~] = PCE_methods.PCE(fc, U, N_ord);

    % define high-fidelity and low-fidelity models
    fcns{1} = @(Z, time_grid) f1(Z, time_grid);   % high-fidelity
    fcns{2} = @(Z) generate_realisations(mean_process, coefficients, basis_index, Z);
    % fcns{3} = @(Z) model3(Z);   % lowest-fidelity

    n_estimate = 300;
    [new_U, ~] = general.generate_legendre_samples(0, d, n_estimate);
    Z = general.generate_inputs(input_ranges, input_means, new_U);
    stats = estimate_statistics(fcns,Z,new_U,vec,time_grid);

    
     G_tot_td = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, true, 'total', dt);
     G_main_td = PCE_methods.calculate_sensitivity_indices(basis_index, coefficients, true, 'main', dt);
    
    % call mfsobol.m with just the high-fidelity model to get Monte
    % Carlo estimate
%     [m,alpha] = optalloc()

    [smmf,stmf,sigsq,~] = mfsobol_vec(fcns(1),d,w(1),stats,m(1),time_grid,vec(1),method);

%     [sigsq', cumsum(sum(  coefficients, 1).^2)];

    a = stats.rho(2)*sqrt(sigsq)./sum(coefficients, 1);

%     mean(a)


%     a = 0.8

    a(a > 1) = 0.95;
    a(a < 0) = 0.05;

    a = 1-a;
    b = 1 - mean(a);

    
    mf_sm(n,:,:) = smmf + a'.*(G_main_td' - smmf);
    mf_st(n,:,:) = stmf + a'.*(G_tot_td' - stmf); 
% 
%     mf_smb(n,:,:) = smmf + b.*(G_main_td' - smmf);

%     mf_sm(n,:,:) = G_main_td' + a*(G_main_td' - smmf);
%     mf_st(n,:,:) = G_tot_td' + a*(G_tot_td' - stmf); 

%     mf_sm(n,:,:)    = sm + alpha(2)*(sm1-sm2);
%     mf_st(n,:,:)    = st + alpha(2)*(st1-st2);
%     % call mfsobol.m with full array of functions to get multifidelity
%     % estimates
%     [sm,st,mu,sigsq] = mfsobol_vec(fcns,d,w,stats,budget,time_grid,vec,method, G_main_td, G_tot_td);
%     avg(n,:,2) = mu; 
%     vr(n,:,2) = sigsq;
%     mf_sm(n,:,:) = sm;
%     mf_st(n,:,:) = st;

    figure(5);
    plot(0:dt:T, reshape(mf_st(n,:,:), N_t,d)); hold on;
    plot(0:dt:T, reshape(mc_st(n,:,:), N_t,d));
    plot(0:dt:T, G_tot_td, 'PC');
    ylim([0,1]); hold off;
n

end

figure(10);
plot(0:dt:T, reshape(mean(mc_sm, 1), N_t, d)', 'PC'); hold on;
plot(0:dt:T, reshape(mean(mf_smb, 1), N_t, d)', 'o');
plot(0:dt:T, reshape(mean(mf_sm, 1), N_t, d)'); hold off;
ylim([0, 1]);


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
title([method,' main sensitivity estimates for Ishigami function'],'interpreter','latex')

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
    title([method,' total sensitivity estimates for Ishigami function'],'interpreter','latex')
end

% % plot variance estimates
% figure(3); clf
% histogram(vr(:,end,1),12,'facecolor',blue); hold on
% histogram(vr(:,end,2),12,'facecolor',red,'facealpha',1)
% legend({'High-fidelity','Multifidelity'},'interpreter','latex'); legend boxoff
% title('Variance estimates for Ishigami function','interpreter','latex')