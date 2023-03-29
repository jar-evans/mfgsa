%% SETUP
clear; close all
addpath('../mfgsa')
addpath('../../td_mfgsa\')

%% set grid
dt = 0.1; T = 250; time_grid = 0:dt:T; N_t = length(time_grid); N_quad = N_t -1;

% uncertain params
bH = 7.5;
kL = 1e6;
X = 168/5;
Z = 70;
g = 7/5;

%% set sample settings
input_mean = [bH, kL, X, Z, g];
input_range = 0.2.*input_mean;

N_p = length(input_mean);
f = @(I) f_cholera(I, 0:dt:T, 2, 45);
% [t, y] = f_cholera(I, 0:dt:T, 2, 45);

d = 5;
% budgets = [100, 250, 750, 1000];

% mse_mct = zeros(4, d); mse_mcm = zeros(4, d);
% mse_pct = zeros(4, d); mse_pcm = zeros(4, d);
% mse_klt = zeros(4, d); mse_klm = zeros(4, d);

% ex_st = [0.260958184278618   0.092823941104023   0.018177267940758   0.097021736546268   0.540445550656443];
% ex_sm = [0.252207054029209   0.088421104950232   0.008822143451783   0.088654590371731   0.536623962706553];
%% COMPUTE MULTIFIDELITY GLOBAL SENSITIVITY ANALYSIS REPLICATES
n_reps = 50;    % number of replicates 

mc_sm = zeros(n_reps,N_quad,d,4);    mc_st = zeros(n_reps,N_quad,d,4);
pc_sm = zeros(n_reps,N_quad,d,4);    pc_st = zeros(n_reps,N_quad,d,4);
kl_sm = zeros(n_reps,d,4);    kl_st = zeros(n_reps,d,4);

method = 'Owen';

% define high-fidelity and low-fidelity models
fcns{1} = @(I, time_grid) f_cholera(I, 0:dt:T, 2, 45);   % high-fidelity
fcns{2} = @(I, time_grid) f_cholera(I, 0:dt:T, 2, 1); % lowest-fidelity
% fcns{3} = @(Z) model3(Z);   % lowest-fidelity
% fcns{1} = @(I) f_cholera(I, time_grid, 2, 45);
w = [1; 0.333];
vec = [1, 3];

% n_estimate = 100;
% new_U = general.generate_legendre_samples(d, n_estimate);
% Z = general.generate_inputs(input_range, input_mean, new_U);
% stats = estimate_statistics(fcns,Z,new_U,vec,0:dt:T);
% 
% for j = 1:length(budgets)
%     p = budgets(j); [m,alpha] = optalloc(p,w,stats)
% end

m = [50, 700];


    for n = 1:n_reps

%         kl_N = 
        

        base = TD_SURROGATE(f, time_grid, input_range, input_mean);
        kl_surrogate = base.generate_KLPC_surrogate(m(1), 15, 3);
        G_tot_kl = kl_surrogate.calculate_sensitivity_indices('total');
        G_main_kl = kl_surrogate.calculate_sensitivity_indices('main');

        kl_I = kl_surrogate.I;
        p = floor(m(2)/(d+2));

        U_A = general.generate_legendre_samples(N_p, p - m(1)); 
        Z_A = general.generate_model_inputs(input_range, input_mean, U_A);

        Z_A = vertcat(kl_I, Z_A);
    
        U_B = general.generate_legendre_samples(N_p, p - m(1)); 
        Z_B = general.generate_model_inputs(input_range, input_mean, U_B);

        Z_B = vertcat(kl_I, Z_B);
    
        yA = fcns{1}(Z_A);
        yB = fcns{1}(Z_B);
    
        yC = zeros(p,N_t,d);
        for i = 1:d
            Z_Ci = Z_B;
            Z_Ci(:,i) = Z_A(:,i);
            yC(:,:,i) = fcns{1}(Z_Ci);
        end
    
        [sm2,st2] = estimate_sobol(method,yA(1:m(1)),yB(1:m(1)),yC(1:m(1)),true, true);
        [sm1,st1] = estimate_sobol(method,yA,yB,yC,true, true);

        alpha = 0.5;
        main = G_main_kl + alpha*(sm1 - sm2)
        tot = G_tot_kl + alpha*(st1 - st2)
    

%         %% establish surrogate
%         N = pc_budgets(j);
% 
% 
%         pc_surrogate = base.generate_PC_surrogate(N, N_ord);
%         G_tot_pc = pc_surrogate.calculate_sensitivity_indices(true, 'total');
%         G_main_pc = pc_surrogate.calculate_sensitivity_indices(true, 'main');
%         G_tot_pc = G_tot_pc';
%         G_main_pc = G_main_pc';

%         mc_sm(n,:,:, j) = sm;
%         mc_st(n,:,:, j) = st;

%         pc_sm(n, :, :, j) = G_main_pc;
%         pc_st(n, :, :, j) = G_tot_pc;

%         mse_pct(j, :) = mse_pct(j, :) + sum(abs(G_tot_pc(end, :) - ex_st)).^2;
%         mse_pcm(j, :) = mse_pcm(j, :) + sum(abs(G_main_pc(end, :) - ex_sm)).^2;


        %% establish surrogate
%         base = TD_SURROGATE(f, time_grid, input_range, input_mean); N_ord = 3; N_KL = 15;
% 
%         N = kl_budgets(j);
% 
%         kl_surrogate = base.generate_KLPC_surrogate(N, N_KL, N_ord);
%         G_tot_kl = kl_surrogate.calculate_sensitivity_indices('total');
%         G_main_kl = kl_surrogate.calculate_sensitivity_indices('main');
%         kl_sm(n, :, j) = G_main_kl; kl_st(n, :, j) = G_tot_kl;
% 
%      
% 
%         mse_klt(j, :) = mse_klt(j, :) + sum(abs(G_tot_kl - ex_st)).^2;
%         mse_klm(j, :) = mse_klm(j, :) + sum(abs(G_main_kl - ex_sm)).^2;
        % error of each s_i is included in same metric
   
    
    end
% end
% 
% mse_mct = (1/n_reps)*mse_mct; mse_mct = mse_mct(:, 1);
% mse_mcm = (1/n_reps)*mse_mcm; mse_mcm = mse_mcm(:, 1);
% 
% mse_pct = (1/n_reps)*mse_pct; mse_pct = mse_pct(:, 1);
% mse_pcm = (1/n_reps)*mse_pcm; mse_pcm = mse_pcm(:, 1);
% 
% mse_klt = (1/n_reps)*mse_klt; mse_klt = mse_klt(:, 1);
% mse_klm = (1/n_reps)*mse_klm; mse_klm = mse_klm(:, 1);
% 
% figure(7); clf
% loglog(mc_budgets, mse_mct', '-o');  hold on;
% loglog(pc_budgets, mse_pct', '-o');
% loglog(kl_budgets, mse_klt', '-o'); grid on;
% legend(flipud(findall(gca,'Tag','Box')), {'MC', 'PCE', 'KLPC'},...
%     'Location','NorthEast','interpreter','latex'); legend boxoff
% bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
% title('Total SI MSE for SIR - 1FHF','interpreter','latex')
% xlabel('Budget (model evaluations)'); ylabel('MSE'); ylim([9e-3, 2]); hold off;
% 
% figure(8); clf
% loglog(mc_budgets, mse_mcm', '-o');  hold on;
% loglog(pc_budgets, mse_pcm', '-o');
% loglog(kl_budgets, mse_klm', '-o'); grid on;
% legend(flipud(findall(gca,'Tag','Box')), {'MC', 'PCE', 'KLPC'},...
%     'Location','NorthEast','interpreter','latex'); legend boxoff
% bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
% title('Main SI MSE for SIR - 1FHF','interpreter','latex')
% xlabel('Budget (model evaluations)'); ylabel('MSE'); ylim([9e-3, 2]); hold off;
% 
% 
% 
% % figure(3);
% % for j = 1:length(budgets)
% %     plot(0:dt:T, reshape(mean(mc_sm(:,:,:,j), 1), N_t, d)'); hold on;  
% % end
% % ylim([0, 1]); hold off;
% 
% 
% %% PLOT ESTIMATOR SPREAD
% warning('off','MATLAB:legend:IgnoringExtraEntries')
% blue = [0       0.4470 0.7410];
% red  = [0.8500  0.3250 0.0908];
% something  = [0.3  0.8 0.2];
% 
% % plot main effect sensitivity indices
% figure(1); clf
% h = boxplot([mc_st(:,end,1,1), mc_st(:,end,1,4),...
%              mc_st(:,end,2,1), mc_st(:,end,2,4), ...
%              mc_st(:,end,3,1), mc_st(:,end,3,4),...
%              mc_st(:,end,4,1), mc_st(:,end,4,4),...
%               mc_st(:,end,5,1), mc_st(:,end,5,4)],...
%     'Colors',[blue; red; blue; red; blue; red;blue; red;blue; red;],'Whisker',10,...
%     'labels',{'$s_t^1$','$s_t^2$','$s_t^3$','$s_t^4$','$s_t^5$',...
%     '$s_t^1$','$s_t^2$','$s_t^3$','$s_t^4$','$s_t^5$'}); hold on;
% set(h,{'linew'},{2}); grid on
% legend(flipud(findall(gca,'Tag','Box')), {'p = 250', 'p = 1000'},...
%     'Location','NorthEast','interpreter','latex'); legend boxoff
% bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
% title('Total SI spread for SIR - 1FHF - MC','interpreter','latex')
% 
% figure(2); clf
% h = boxplot([mc_sm(:,end,1,1), mc_sm(:,end,1,4),...
%              mc_sm(:,end,2,1), mc_sm(:,end,2,4), ...
%              mc_sm(:,end,3,1), mc_sm(:,end,3,4),...
%              mc_sm(:,end,4,1), mc_sm(:,end,4,4),...
%               mc_sm(:,end,5,1), mc_sm(:,end,5,4)],...
%     'Colors',[blue; red; blue; red; blue; red;blue; red;blue; red;],'Whisker',10,...
%     'labels',{'$s_m^1$','$s_m^2$','$s_m^3$','$s_m^4$','$s_m^5$',...
%     '$s_m^1$','$s_m^2$','$s_m^3$','$s_m^4$','$s_m^5$'}); hold on;
% set(h,{'linew'},{2}); grid on
% legend(flipud(findall(gca,'Tag','Box')), {'p = 250', 'p = 1000'},...
%     'Location','NorthEast','interpreter','latex'); legend boxoff
% bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
% title('Main SI spread for SIR - 1FHF - MC','interpreter','latex')
% 
% 
% 
% 
% 
% 
