function [time_grid, Y] = f_cholera(I, dt, T, ic, tols, N, QoI)



    size(I)
    time_grid = 0:dt:T;

    tspan = [0 T];%[5e-2, 2e2];    
    opts = odeset('RelTol', tols(2), 'AbsTol', tols(1));
    for i = 1:N

%         if N == 1
%             [t,y] = ode45(@(t,y) cholera(t, y, ic(1), I), tspan, ic(2:length(ic)), opts);
%         else
            [t,y] = ode45(@(t,y) cholera(t, y, ic(1), I(i, :)), tspan, ic(2:length(ic)), opts);
%         end

        y = y(:, QoI)';

        if i == 1
            Y = fix_dt(t, y, time_grid);
        else
%             size(Y)
%             size(y)
%             size(fix_dt(t, y, time_grid))
            Y = vertcat(Y, fix_dt(t, y, time_grid));
        end

    end

end

function y = fix_dt(t, y, time_grid)
    r = interp1(time_grid, time_grid, t, 'nearest');
    [~, filt] = unique(r);
    y = y(filt);
end

function dydt = cholera(t, y, N_pop, p)
    dydt = zeros(5,1);
%     y(1) = S 
%     y(2) = I
%     y(3) = R
%     y(4) = B_H
%     y(5) = B_L
    dydt(1) = p(5)*N_pop - p(1)*y(1)*y(5)/(p(3) + y(5)) - p(2)*y(1)*y(4)/(p(4) + y(4)) - p(5)*y(1);
    dydt(2) = p(1)*y(1)*y(5)/(p(3) + y(5)) + p(2)*y(1)*y(4)/(p(4) + y(4)) - (p(9) + p(5))*y(2);
    dydt(3) = p(9)*y(2) - p(5)*y(3);
    dydt(4) = p(7)*y(2) - p(6)*y(4);
    dydt(5) = p(6)*y(4) - p(8)*y(5);
end


