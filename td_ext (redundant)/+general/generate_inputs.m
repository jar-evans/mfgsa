function I = generate_inputs(ranges, means, legendre_samples)

    [N, N_p] = size(legendre_samples);
    assert(numel(ranges) == numel(means) && numel(means) == N_p);

    I = zeros(N, N_p);

    for i = 1:N_p
        I(1:N, i) = 0.5*ranges(i)*legendre_samples(1:N, i) + means(i);
    end
end

% A_range = 0.25;
% B_range = 1.25;
% L_range = 0.5;
% 
% A = (0.125*U(1:N, 1)) + 0.5;
% B = ((5/8)*U(1:N, 2)) + (25/8);
% L = (0.25*U(1:N, 3)) - 1;
% 
% Af = 0.5*ones(N,1);
% Bf = (25/8)*ones(N,1);
% Lf = -1*ones(N,1);
% U_af = 0; U_bf = 0; U_lf = 0;