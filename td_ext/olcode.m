% figure(6);
% plot(1-p, 'PC'); hold on;
% plot(y.*y);
% plot(tot); hold off;
% 
% G = 1 - (sum(p)/sum(tot))





% Ki_filt
% 
% fuck_it = sum(proj_matrix .* proj_matrix, 2);
% 
% c_p = c(Ki_filt, :);
% fuck_it_p = fuck_it(Ki_filt, :);
% 
% w = 1;
% c2_p = w .* c_p .* c_p;
% c2_tot = w .* c .* c;
% 
% a = sum(c2_p, 1)*dt;
% b = sum(c2_tot, 1)*dt;
% 
% 
% size(a)
% R = 20;
% BB = zeros(1, R);
% AA = zeros(1, R);
% 
% for r = R:-1:1
%     AA(r) = sum(a(1:int64(length(a)/r)));
%     BB(r) = sum(b(1:int64(length(b)/r)));
% end
% 
% plot(1:R, AA, 1:R, BB);





% num = c2_p .* fuck_it_p;
% den = c2_tot .* fuck_it;
% 
% G = sum(sum(num, 2))/sum(sum(den, 2))

% plot(1:N_quad, sum(y))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = 1;%N_p/N;%1e-10;
% K_B = basis_index(basis_index(:, 1) > 0, :);
% 
% 
% sum_p = 0;
% for n = 1:length(K_B)
%     for m = 1:N_quad
%         temp = 0;
%         for k = 1:N_p
%             res = legendre(K_B(n, k), U(:, k)).^2;
%             if (numel(res) > 1)
%                 sum_p = sum_p + (w * c(:, m) .* c(:, m) * res(1, :));
%             else
%                 sum_p = sum_p + (w * c(:, m) .* c(:, m) * res);
%             end
%         end
%         sum_p = sum_p + temp;
% 
%     end
% end

% sub_sums = zeros(1, N_PC);
% 
% sum_t = 0;
% for n = 1:N_PC
%     temp2 = 0;
%     for m = 1:N_quad
%         temp = 1;
%         for k = 1:N_p
%             size(temp);
%             res = legendre(basis_index(n, k), U(:, k), 'norm');
%             if (numel(res) > 1)
%                 temp = temp .* res(1, :) .* res(1, :);
% %                 temp = temp + (w .* c(:, m) .* c(:, m) .* res(1, :));
%             else
%                 temp = temp .* res .* res;
%            
% %                 temp = temp + (w .* c(:, m) .* c(:, m) * res);
%             end
% 
%             
% 
%         end
%         
%         temp2 = temp2 + (temp .* c(:, m) .* c(:, m));
%     end
%     size(temp2);
%     sub_sums(n) = temp2;
% end
% 
% sub_sums;


% sum(sum_p)/sum(sum_t)
% plot(1:35, sum_p); hold on; plot(1:35, sum_t); hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res = 0; temp = zeros(1, N_quad); temp2 = 0; temp3 = 0;

K_L = basis_index(basis_index(:, 3) > 0, :);
% for n = 1:length(K_L)
for n = 1:N_PC
    temp = c(n, :);% .* c(n, :);
    for k = 1:N_p        
        norm = 2/((2*basis_index(n, k)) + 1);
        temp2 = legendre(basis_index(n, k), U(:, k))/norm; 
%         temp2 = legendre(K_L(n, k), leg_samples(:, k), 'norm'); 
        temp3 = temp2(1, :);% .* temp2(1, :);        
    end
    temp = transpose(temp) * temp3;
    res = res + temp;
end

figure(11);
res = transpose(res);
plot(sum(res)); hold on;
plot(1:N_quad, res)
plot(res(1,:));

figure(12);
plot(max(res)); hold on; plot(min(res)); plot(mean(res)); hold off;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% alg 3 stuff down there


% L_sort = sort(U_lf);
% KUc = [U_lf, U_bf];
% 
% res = zeros(I, N_PC_quad);
% 
% for i = 1:I
%     for n = 1:N_PC
% %         res(i, :) = res(i, :) + c(i, n);
%         temp2 = c(i, n);
%         for k = 1:N_p-1
%             temp = legendre(basis_index(n, k), KUc(k));
% %             res(i, :) = res(i, :) .* temp(1);
%             temp2 = temp2 .* temp(1);
%         end
%         if (basis_index(n, k) == 0)
%             temp = 0;
%         else
%             temp = legendre(basis_index(n, k), L_sort);
%         end
% %         res(i, :) = res(i, :) .* temp(1, :);
%         temp2 = temp2 .* temp(1,:);
%         res(i, :) = res(i, :) + temp2;
%     end
% end



res = zeros(I, N_PC_quad);

for i = 1:I
    for n = 1:N_PC
%         res(i, :) = res(i, :) + c(i, n);
        temp2 = c(i, n);
        for k = 1:N_p
            temp = legendre(basis_index(n, k), U(k), 'norm');
%             res(i, :) = res(i, :) .* temp(1);
            temp2 = temp2 .* temp(1);
        end

%         res(i, :) = res(i, :) .* temp(1, :);
        res(i, :) = res(i, :) + temp2;
    end
end



fuck = transpose(res)*transpose(e);
figure(7); hold on;
plot(transpose(sum(fc))); plot(sum(fuck), "PC"); hold off;