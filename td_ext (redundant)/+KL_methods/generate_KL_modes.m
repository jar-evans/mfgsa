function fi = generate_KL_modes(fc, e)

    N = size(fc, 1);
    N_KL = size(e, 2);
    fi = zeros(N_KL, N);

    w = 1/N;%ones(1, N_quad);%double(1/N_KL); %
    
    for i = 1:N_KL
        fi(i, :) = w*(fc*e(:,i)); % is the use of w = dt right here?   % this only works for uniform quadrature
    end

end

