function [sm,st] = estimate_sobol(method,yA,yB,X,vec)%,gen)
% INPUTS
% method    'Owen' or 'Saltelli' or 'Gamboa'
% yA        N-by-1 first set of function evaluations
% yB        for 'Owen' or 'Saltelli' only: N-by-1 second set of function 
%           evaluations. For 'Gamboa' pass an empty array []
% X         **for 'Owen' or Saltelli': N-by-d function evaluationsthe i-th 
%           column of X was evaluated at the same inputs that led to yB 
%           except replacing the i-th column of the input with the i-th 
%           column of yA. 
%           **for 'Gamboa': N-by-d array of inputs that led to outputs yA
%
% OUTPUTS
% sm        d-by-1 vector of main effect sensitivities
% st        d-by-1 vector of total effect sensitivities (zero output for
%           method = 'Gamboa')
%
% AUTHORS
% Elizabeth Qian (www.elizabethqian.com)
% Giuseppe Cataldo
%
% LAST UPDATED
% 18 May 2022

[N,d] = size(X);

muA  = mean(yA);
varA = var(yA);

if ~vec
    switch method
        case 'Owen'
            sm = 2*N/(2*N-1)*(yA'*X/N - (muA+mean(X)).^2/4 + (varA+var(X))/(4*N))/varA;
            st = 1/(2*N)*sum((yB-X).^2)/varA;
        case 'Saltelli'
            sm = (1/(N-1)*yA'*X - muA^2)/varA;
            st = 1 - (1/(N-1)*yB'*X - muA^2)/varA;
        case 'Gamboa'
            % Get X ordering
            [~,px] = sort(X);
            % Get pi_j
            [~,pi_j] = sort(px);
            % Get N_j with shift
            argpiinv = mod(pi_j,N)+1;
            % Index in x of thing with this rank
            N_j = zeros(N,d);
            for i = 1:d
                N_j(:,i) = px(argpiinv(:,i),i);
            end
            % Get output corresponding to new rank
            y_Nj = yA(N_j);
            sm = (1/N*yA'*y_Nj - muA^2)/varA;
            st = zeros(1,d);
    end
else
    
    den = varA;
    
    %% maybe its worth having a pw option s.t. pw variance of the process can
    % be visualised and compared back to the td paper
    switch method
        case 'Owen'
            [N,N_t,d] = size(X);
    
            a = permute(X, [1 3 2]);        
            resA = matmul_2by3(yA', a, N_t, d);
    
            meanX = permute(mean(X,1), [2 3 1]);
            varX = permute(var(X,1), [2 3 1]);
    
            num_m = 2*N/(2*N-1)*(resA/N - (muA'+meanX).^2/4 + (varA'+varX)/(4*N));
            num_t = 1/(2*N)*sum((yB-X).^2);

            den = reshape(den, N_t, 1);
            num_m = reshape(num_m, N_t, d);
            num_t = reshape(num_t, N_t, d);

%             if gen

                num_m_quad = diff(num_m,1)/2 + num_m(1:end-1,:);
                num_t_quad = diff(num_t,1)/2 + num_t(1:end-1,:);
                den_quad = diff(den,1)/2 + den(1:end-1,:);

                sm = cumsum(num_m_quad)./cumsum(den_quad);
                st = cumsum(num_t_quad)./cumsum(den_quad);

%             else
%                 sm = num_m./den;
%                 st = num_t./den;
% 
%             end
    
        case 'Saltelli' 
            [N,N_t,d] = size(X);
    
            a = permute(X, [1 3 2]);
            
            resA = matmul_2by3(yA', a, N_t, d);
            resB = matmul_2by3(yB', a, N_t, d);
    
            num_m = (1/(N-1)*resA - muA'.^2);
            num_t = (1/(N-1)*resB - muA'.^2);
    
            sm = cumsum(num_m)./cumsum(den)';
            st = 1 - (cumsum(num_t)./cumsum(den)');
    
        case 'Gamboa'
            error("NotImplementedError");
    end
end
end

function res = matmul_2by3(b, a, I, J)

    res = zeros(I,J);
    for i = 1:I
        for j = 1:J
            res(i, j) = b(i,:) * a(:,j,i);
        end
    end
end

