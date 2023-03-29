function y = generate_approximation(projection_matrix, coefficients, basis_filter)
    
    if ~exist('basis_filter', 'var')
        y = sum(transpose(projection_matrix)*coefficients);
    else
        y = sum(transpose(projection_matrix(basis_filter, :))* ...
            coefficients(basis_filter, :)); 
    end
    
end

%% is this right though?
%    is the approx just c*proj, proj includes a quadrature weight/norm ...
%    should it not just be c*proj*norm/v ???
