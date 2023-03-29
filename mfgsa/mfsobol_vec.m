function [sm,st,mu,sigsq] = mfsobol_vec(fcns,d,w,stats,p,time_grid,vec,estim, input_mean,input_range)
% computes multifidelity estimate of mean, variance, and Sobol' main and
% total effect sensitivity indices
%
% INPUTS
% fcns      k-by-1 cell array of functions corresponding to different models
% d         dimension of uncertaint input
% w         k-by-1 vector of computational costs for functions in fcns
% stats     struct containing statistics of models in fcns
% p         total computational budget
% vec       k-by-1 vector that indicates whether the models in fcns are
%           vectorized or not. 0 (false) is default. 1 is vectorized. all
%           2's indicates that functions bootstrap from samples.
% estim     'Owen' or 'Saltelli' or 'Gamboa' argument to be passed to 
%           estimate_sobol
%
% OUTPUTS
% mu        multifidelity mean estimate of high-fidelity model in fcns{1}
% sigsq     multifid. variance estimate of high-fidelity model in fcns{1}
% sm        d-by-1 vector of Sobol' main effect sensitivity indices
% st        d-by-1 vector of Sobol' total effect sensitivity indices
%
% AUTHORS
% Elizabeth Qian (elizqian@alum.mit.edu) 
% Giuseppe Cataldo 
%
% LAST UPDATED
% 18 May 2022

N_t = length(time_grid);

if nargin <= 5
    vec = zeros(size(fcns));
end

if any(vec == 2)
    assert(all(vec==2),'If computing indices via bootstrapping, fcns need to have special form and the vectorization flag needs to be 2 for all models')
end

fixing = ~strcmp(estim,'Gamboa'); % flag for whether method needs fixing inputs

% get optimal number of evaluations and weights 
if fixing % use effective budget
    [m,alpha] = optalloc(p/(d+2),w,stats);
    % optalloc(p,w,stats); m(1) = floor(m(1)/(d+2));
else % use whole budget
    [m,alpha] = optalloc(p,w,stats);
end

% if length(vec) ==1 
%     m = 50;
%     alpha = 1;
% else
% m = [15 35];
% alpha = [1 0.5];
% end

% get two sets of independent inputs (generate_inputs.m for problem must be
% on MATLAB search path)
if vec(1) ~= 2 
    % if not bootstrapping, generate new samples
    U_A = general.generate_legendre_samples(d, m(end));
    Z_A = general.generate_inputs(input_range, input_mean, U_A);

    if fixing
        U_B = general.generate_legendre_samples(d, m(end));
        Z_B = general.generate_inputs(input_range, input_mean, U_A);
    end
else
    % if bootstrapping, just pull a random set of indices
    error("Bootstrapping - NotImplementedError")
    Z_ind = generate_inputs(m(end));
end



if vec(1) == 2
    % if bootstrapping, pull all outputs from high-fi samples now
    error("Bootstrapping - NotImplementedError")
    [yA, yB, yC] = fcns{1}(Z_ind(1:m(1)));
end

% compute all evaluations of high-fidelity model
if vec(1) ~= 0
    % if model is vectorized, evaluate at all inputs at once
    yA = fcns{1}(Z_A(1:m(1),:), time_grid);
    if fixing
        yB = fcns{1}(Z_B(1:m(1),:), time_grid);
    end
elseif vec(1) == 0 
    % if model is not vectorized, evaluate inputs in loop
    yA = zeros(m(1),N_t); 
    for j = 1:m(1)
        yA(j,:) = fcns{1}(Z_A(j,:), time_grid);
        
    end

    if fixing
        yB = zeros(m(1),N_t);
        for j = 1:m(1)
            yB(j,:) = fcns{1}(Z_B(j,:), time_grid);
        end
    end
end

if fixing 
    if vec(1) ~= 2
        % if not bootstrapping
        yC = zeros(m(1),N_t,d);
        for i = 1:d
            Z_Ci = Z_B(1:m(1),:);
            Z_Ci(:,i) = Z_A(1:m(1),i);
    
            if vec(1)
                yC(:,:,i) = fcns{1}(Z_Ci, time_grid);
            else
                for j = 1:m(1)
                    yC(j,:,i) = fcns{1}(Z_Ci(j,:), time_grid);
                end
            end
        end
    end
end

% initialize all statistics with their high-fidelity values
if fixing
    mu      = mean([yA; yB]); assert(length(mu) == N_t);
    sigsq   = var([yA; yB]);
    [sm,st] = estimate_sobol(estim,yA,yB,yC,true);
else
    mu      = mean(yA);
    sigsq   = var(yA);
    [sm,st] = estimate_sobol(estim,yA,[],Z_A(1:m(1),:),true);
end

% loop through low-fidelity models
for j = 2:length(m)

    if vec(j) ~= 3
    
        if vec(j) == 2
            error("Bootstrapping - NotImplementedError")
            [yA, yB, yC] = fcns{j}(Z_ind(1:m(j)));
        else
            % get function evalutions
            if vec(j) == 1
                yA = fcns{j}(Z_A(1:m(j),:), time_grid);
                if fixing
                    yB = fcns{j}(Z_B(1:m(j),:), time_grid);
                end
            elseif vec(j) == 3
                
                yA = fcns{j}(U_A(1:m(j),:));
                if fixing
                    yB = fcns{j}(U_B(1:m(j),:));
                end
            else
                yA = zeros(m(j),N_t,1); 
                for k = 1:m(j)
                    yA(k,:) = fcns{j}(Z_A(k,:), time_grid);
                end
                if fixing
                    yB = zeros(m(j),N_t,1);
                    for k = 1:m(j)
                        yB(k,:) = fcns{j}(Z_B(k,:), time_grid);
                    end
                end
            end
            
            if fixing
                yC = zeros(m(j),N_t,d);
                for i = 1:d
                    Z_Ci = Z_B(1:m(j),:);
                    Z_Ci(:,i) = Z_A(1:m(j),i);
        
                    if vec(j) == 1
                  
                        yC(:,:,i) = fcns{j}(Z_Ci, time_grid);
                    elseif vec(j) == 0
                        for k = 1:m(j)
                            yC(k,:,i) = fcns{j}(Z_Ci(k,:), time_grid);
                        end
                    elseif vec(j) == 3
                        U_Ci = U_B(1:m(j),:);
                        U_Ci(:,i) = U_A(1:m(j),i);
                  
                        yC(:,:,i) = fcns{j}(U_Ci);
    
                    end
                end
            end
        end
    end
   
    if vec(j) == 3

        if fixing
            mu    = 0;
            sigsq = 0;
            [sm1,st1] = estimate_sobol(estim,yA(1:m(j-1),:),yB(1:m(j-1),:), yC(1:m(j-1),:,:),true);
            [sm2,st2] = estimate_sobol(estim,yA(1:m(j-1),:),yB(1:m(j-1),:), yC(1:m(j-1),:,:),true);
        else
            mu    = 0;
            sigsq = 0;
            [sm1,st1] = estimate_sobol(estim,yA(1:m(j),:),[],Z_A(1:m(j),:),true);
            [sm2,st2] = estimate_sobol(estim,yA(1:m(j-1),:),[], Z_A(1:m(j-1),:),true);
        end

        sm    = sm + alpha(j)*(sm1-sm2);
        st    = st + alpha(j)*(st1-st2);


    else

    % add low-fi correction to existing estimate
        if fixing
            
            mu    = mu + alpha(j)*(mean([yA; yB]) - mean([yA(1:m(j-1),:); yB(1:m(j-1),:)]));
            sigsq = sigsq + alpha(j)^2*(var([yA; yB]) - var([yA(1:m(j-1),:); yB(1:m(j-1),:)]));
            [sm1,st1] = estimate_sobol(estim,yA,yB,yC,true);
            [sm2,st2] = estimate_sobol(estim,yA(1:m(j-1),:),yB(1:m(j-1),:), yC(1:m(j-1),:,:),true);
        else
            mu    = mu + alpha(j)*(mean(yA) - mean(yA(1:m(j-1),:)));
            sigsq = sigsq + alpha(j)^2*(var(yA) - var(yA(1:m(j-1),:)));
            [sm1,st1] = estimate_sobol(estim,yA(1:m(j),:),[],Z_A(1:m(j),:),true);
            [sm2,st2] = estimate_sobol(estim,yA(1:m(j-1),:),[], Z_A(1:m(j-1),:),true);
        end
        sm    = sm + alpha(j)*(sm1-sm2);
        st    = st + alpha(j)*(st1-st2);

    end
end