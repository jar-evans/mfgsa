%% SETUP
clear; close all
addpath('../mfgsa')
addpath('../../td_mfgsa\')

%% generate realisations
dt = 0.1; T = 10; time_grid = 0:dt:T; N_t = length(time_grid);
input_range = [0.25 1.25 0.5]; input_mean = [0.5 25/8 -1];


d = 3; budgets = [100, 1000, 10000];

mse_mct = zeros(length(budgets), d); mse_mcm = zeros(length(budgets), d);
mse_pct = zeros(length(budgets), d); mse_pcm = zeros(length(budgets), d);
mse_klt = zeros(length(budgets), d); mse_klm = zeros(length(budgets), d);

ex_st = [0.052116592012602   0.871787056655515   0.127187780429433];
ex_sm = [0.019649685321925   0.822243265739821   0.109172883919503];
%% COMPUTE MULTIFIDELITY GLOBAL SENSITIVITY ANALYSIS REPLICATES
n_reps = 100;    % number of replicates 

mc_sm = zeros(n_reps,N_t-1,d,length(budgets));    mc_st = zeros(n_reps,N_t-1,d,length(budgets));
pc_sm = zeros(n_reps,N_t-1,d,length(budgets));    pc_st = zeros(n_reps,N_t-1,d,length(budgets));
kl_sm = zeros(n_reps,d,length(budgets));    kl_st = zeros(n_reps,d,length(budgets));

method = 'Owen';

fcns{1} = @(Z) f2(Z, time_grid);   % high-fidelity
f = @(I) f2(I, time_grid);

for j = 1:length(budgets)

    p = budgets(j)/(d+2);
    for n = 1:n_reps

        n
    
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
    
        [sm,st] = estimate_sobol(method,yA,yB,yC,true);
    
        mc_sm(n,:,:, j) = sm;
        mc_st(n,:,:, j) = st;

        %% establish surrogate
        base = TD_SURROGATE(f, time_grid, input_range, input_mean); N_ord = 4; N_KL = 8;

        N = budgets(j);

        pc_surrogate = base.generate_PC_surrogate(N, N_ord);
        G_tot_pc = pc_surrogate.calculate_sensitivity_indices(true, 'total');
        G_main_pc = pc_surrogate.calculate_sensitivity_indices(true, 'main');
        G_tot_pc = G_tot_pc';
        G_main_pc = G_main_pc';
        pc_sm(n, :, :, j) = G_main_pc; pc_st(n, :, :, j) = G_tot_pc;

        kl_surrogate = base.generate_KLPC_surrogate(N, N_KL, N_ord);
        G_tot_kl = kl_surrogate.calculate_sensitivity_indices('total');
        G_main_kl = kl_surrogate.calculate_sensitivity_indices('main');
        kl_sm(n, :, j) = G_main_kl; kl_st(n, :, j) = G_tot_kl;

        mse_mct(j, :) = mse_mct(j, :) + sum(abs(st(end, :) - ex_st)).^2;
        mse_mcm(j, :) = mse_mcm(j, :) + sum(abs(sm(end, :) - ex_sm)).^2;

        mse_pct(j, :) = mse_pct(j, :) + sum(abs(G_tot_pc(end, :) - ex_st)).^2;
        mse_pcm(j, :) = mse_pcm(j, :) + sum(abs(G_main_pc(end, :) - ex_sm)).^2;

        mse_klt(j, :) = mse_klt(j, :) + sum(abs(G_tot_kl - ex_st)).^2;
        mse_klm(j, :) = mse_klm(j, :) + sum(abs(G_main_kl - ex_sm)).^2;
        % error of each s_i is included in same metric
   
    
    end
end

mse_mct = (1/n_reps)*mse_mct; mse_mct = mse_mct(:, 1);
mse_mcm = (1/n_reps)*mse_mcm; mse_mcm = mse_mcm(:, 1);

mse_pct = (1/n_reps)*mse_pct; mse_pct = mse_pct(:, 1);
mse_pcm = (1/n_reps)*mse_pcm; mse_pcm = mse_pcm(:, 1);

mse_klt = (1/n_reps)*mse_klt; mse_klt = mse_klt(:, 1);
mse_klm = (1/n_reps)*mse_klm; mse_klm = mse_klm(:, 1);

figure(4); clf
loglog(budgets, mse_mct', '-o');  hold on;
loglog(budgets, mse_pct', '-o');
loglog(budgets, mse_klt', '-o'); grid on;
legend(flipud(findall(gca,'Tag','Box')), {'MC', 'PCE', 'KLPC'},...
    'Location','NorthEast','interpreter','latex'); legend boxoff
bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
title('Total SI MSE for mech. oscillator - 1FLF','interpreter','latex')
xlabel('Budget (model evaluations)'); ylabel('MSE'); hold off;

figure(5); clf
loglog(budgets, mse_mcm', '-o');  hold on;
loglog(budgets, mse_pcm', '-o');
loglog(budgets, mse_klm', '-o'); grid on;
legend(flipud(findall(gca,'Tag','Box')), {'MC', 'PCE', 'KLPC'},...
    'Location','NorthEast','interpreter','latex'); legend boxoff
bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
title('Main SI MSE for mech. oscillator - 1FLF','interpreter','latex')
xlabel('Budget (model evaluations)'); ylabel('MSE'); hold off;



% figure(3);
% for j = 1:length(budgets)
%     plot(0:dt:T, reshape(mean(mc_sm(:,:,:,j), 1), N_t, d)'); hold on;  
% end
% ylim([0, 1]); hold off;


%% PLOT ESTIMATOR SPREAD
warning('off','MATLAB:legend:IgnoringExtraEntries')
blue = [0       0.4470 0.7410];
red  = [0.8500  0.3250 0.0908];
something  = [0.3  0.8 0.2];

% plot main effect sensitivity indices
figure(1); clf
h = boxplot([mc_st(:,end,1,1), mc_st(:,end,1,2), mc_st(:,end,1,3),...
             mc_st(:,end,2,1), mc_st(:,end,2,2), mc_st(:,end,2,3), ...
             mc_st(:,end,3,1), mc_st(:,end,3,2), mc_st(:,end,3,3)],...
    'Colors',[blue; red; something; blue; red; something; blue; red; something],'Whisker',10,...
    'labels',{'$s_t^1$','$s_t^2$','$s_t^3$', '$s_t^1$','$s_t^2$','$s_t^3$', '$s_t^1$','$s_t^2$','$s_t^3$'}); hold on;
set(h,{'linew'},{2}); grid on
legend(flipud(findall(gca,'Tag','Box')), {'p = 100', 'p = 1000', 'p = 10000'},...
    'Location','NorthEast','interpreter','latex'); legend boxoff
bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
title('Total SI spread for mech. oscillator - 1FLF - MC','interpreter','latex')


figure(2); clf
h = boxplot([pc_st(:,end,1,1), pc_st(:,end,1,2), pc_st(:,end,1,3),...
             pc_st(:,end,2,1), pc_st(:,end,2,2), pc_st(:,end,2,3), ...
             pc_st(:,end,3,1), pc_st(:,end,3,2), pc_st(:,end,3,3)],...
    'Colors',[blue; red; something; blue; red; something; blue; red; something],'Whisker',10,...
    'labels',{'$s_t^1$','$s_t^2$','$s_t^3$', '$s_t^1$','$s_t^2$','$s_t^3$', '$s_t^1$','$s_t^2$','$s_t^3$'}); hold on;
set(h,{'linew'},{2}); grid on
legend(flipud(findall(gca,'Tag','Box')), {'p = 100', 'p = 1000', 'p = 10000'},...
    'Location','NorthEast','interpreter','latex'); legend boxoff
bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
title('Total SI spread for mech. oscillator - 1FLF - PCE','interpreter','latex')

figure(3); clf
h = boxplot([kl_st(:,1,1), kl_st(:,1,2), kl_st(:,1,3),...
             kl_st(:,2,1), kl_st(:,2,2), kl_st(:,2,3), ...
             kl_st(:,3,1), kl_st(:,3,2), kl_st(:,3,3)],...
    'Colors',[blue; red; something; blue; red; something; blue; red; something],'Whisker',10,...
    'labels',{'$s_m^1$','$s_m^2$','$s_m^3$', '$s_m^1$','$s_m^2$','$s_m^3$', '$s_m^1$','$s_m^2$','$s_m^3$'}); hold on;
set(h,{'linew'},{2}); grid on
legend(flipud(findall(gca,'Tag','Box')), {'p = 100', 'p = 1000', 'p = 10000'},...
    'Location','NorthEast','interpreter','latex'); legend boxoff
bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
title('Total SI spread for mech. oscillator - 1FLF - KLPC','interpreter','latex')
% end



