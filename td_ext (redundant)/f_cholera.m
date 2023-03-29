function [Y, time_grid] = f_cholera(I, time_grid, QoI, solver)

    N_pop = 10000;
    S0 = N_pop - 1; I0 = 1; R0 = 0;
    B_H0 = 0; B_L0 = 0;
    tol_abs = 1e-6;
    tol_rel = 1e-6;

    N = size(I, 1);
    
    
    ic = [N_pop, S0, I0, R0, B_H0, B_L0];
    tols = [tol_abs, tol_rel];

    q = [1.5, 1/1560, 7/30];

    switch solver
        case 45

%             if N < 100
%                 test_i = 1:N;
%             else
%                 test_i = 1:100;
%             end
                
            tspan = [min(time_grid), max(time_grid)];%[5e-2, 2e2];
            opts = odeset('RelTol', tols(2), 'AbsTol', tols(1));    

%             met = zeros(1,100);
%             for i = test_i
%                 [t,~] = ode45(@(t,y) cholera(t, y, ic(1), I(i, :), q), tspan, ic(2:length(ic)), opts);
%                 met(i) = length(t);
%             end
% 
%             grid = min(met);
%             a = floor(grid/100);
%             time_grid = linspace(0, 250, 100*a); length(time_grid)
        
            Y = zeros(N, length(time_grid));
            for i = 1:N
        
                [t,y] = ode45(@(t,y) cholera(t, y, ic(1), I(i, :), q), tspan, ic(2:end), opts);
  
                y = y(:, QoI)';
                Y(i, :) = fix_dt(t, y, time_grid);
            end
    
        case 4
            Y = zeros(N, length(time_grid));
            for i = 1:N
                y = ode4(@(t,y) cholera(t, y, ic(1), I(i, :), q), time_grid, ic(2:end));
                y = y(:, QoI)';
                Y(i, :) = fix_dt(t, y, time_grid);%y;
            end
        case 1
            Y = zeros(N, length(time_grid));
            for i = 1:N
                y = ode1(@(t,y) cholera(t, y, ic(1), I(i, :), q), 0:0.05:250, ic(2:end));
                y = y(:, QoI)';
                Y(i, :) = fix_dt(0:0.05:250, y, time_grid);%y;
            end
    
    end
end
%% this function is not to lower the fidelity -> to ensure the realisations 
%  have consistent lengths, s.t. future calculations are not over complicated
function y = fix_dt(t, y, time_grid) 
%     r = interp1(t, t, time_grid, 'nearest');
%     [~, filt] = unique(r);
%     y = y(filt);

    y = interp1(t, y, time_grid);
end

function dydt = cholera(t, y, N_pop, p, q)

% p - uncertain params [bH, kL, x, z, g]
% q - fixed params [bL, b, d]
    kH = p(2)/700;
    dydt = zeros(5,1);
%     y(1) = S 
%     y(2) = I
%     y(3) = R
%     y(4) = B_H
%     y(5) = B_L
    dydt(1) = q(2)*N_pop - q(1)*y(1)*y(5)/(p(2) + y(5)) - p(1)*y(1)*y(4)/(kH + y(4)) - q(2)*y(1);
    dydt(2) = q(1)*y(1)*y(5)/(p(2) + y(5)) + ...
        p(1)*y(1)*y(4)/(kH + y(4)) - ...
        (p(5) + q(2))*y(2);
    dydt(3) = p(5)*y(2) - q(2)*y(3);
    dydt(4) = p(4)*y(2) - p(3)*y(4);
    dydt(5) = p(3)*y(4) - q(3)*y(5);
end

function Y = ode4(odefun,tspan,y0,varargin)
%ODE4  Solve differential equations with a non-adaptive method of order 4.
%   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the classical Runge-Kutta method of order 4.   
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode4(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%

if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

try
  f0 = feval(odefun,tspan(1),y0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
  error(msg);  
end  

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end  

neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);
F = zeros(neq,4);

Y(:,1) = y0;
for i = 2:N
  ti = tspan(i-1);
  hi = h(i-1);
  yi = Y(:,i-1);
  F(:,1) = feval(odefun,ti,yi,varargin{:});
  F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),varargin{:});
  F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),varargin{:});  
  F(:,4) = feval(odefun,tspan(i),yi+hi*F(:,3),varargin{:});
  Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
end
Y = Y.';
end

function Y = ode1(odefun,tspan,y0,varargin)
%ODE1  Solve differential equations with a non-adaptive method of order 1.
%   Y = ODE1(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE1(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN.
%   The solver implements the forward Euler method of order 1.   
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode1(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%

if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

try
  f0 = feval(odefun,tspan(1),y0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
  error(msg);  
end  

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end  

neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);

Y(:,1) = y0;
for i = 1:N-1 
  Y(:,i+1) = Y(:,i) + h(i)*feval(odefun,tspan(i),Y(:,i),varargin{:});
end
Y = Y.';
end





