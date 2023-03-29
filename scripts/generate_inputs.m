function Z = generate_inputs(N)
    U = generate_legendre_samples(3, N); 
    A = (0.125*U(1:N, 1)) + 0.5;
    B = ((5/8)*U(1:N, 2)) + (25/8);
    L = (0.25*U(1:N, 3)) - 1;

    Z = horzcat(A, B, L);
end


function U = generate_legendre_samples(N_p, N)
        U = 2*rand(N, N_p) - 1;
end