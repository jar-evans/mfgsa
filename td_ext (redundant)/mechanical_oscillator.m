import PCE_methods.*; import general.*;

T = 10;
dt = .01;
N_quad = int64(T/dt) + 1; %length(0:dt:T); % NOTE: should i change N_quad to
                                           % number of intervals and then
                                           % introduce new var N_points

%% generate realisations
input_ranges = [0.25 1.25 0.5];
input_means = [0.5 25/8 -1];

N_p = 3; grid_h = 0; N = 1000;

% [U, N] = general.generate_legendre_samples(grid_h, N_p);
[U, N] = general.generate_legendre_samples(grid_h, N_p, N); 
I = general.generate_inputs(input_ranges, input_means, U);

%% generate N realisations of f, and then centre f
[fmk, time_grid] = f_mechosc(I, T, dt);

%% set time grid % time discr may not always be unique
N_t = length(time_grid);    % #points in time discretisation
N_h = N_t - 1;               % #intervals in time discretisation

G_alg2 = gen_alg2(fmk, time_grid, U, 4, true, 'total');

% G_alg3 = gen_alg3(fmk, time_grid, U, 1e-5, 4, 'main');

