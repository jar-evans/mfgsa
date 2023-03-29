function c = calculate_PCE_coefficients(projection_matrix, function_evaluations)

    N_PC = size(projection_matrix, 1);
    N_quad = size(function_evaluations, 2);

    c = zeros(N_PC, N_quad);
    for tm = 1:N_quad
        c(:, tm) = projection_matrix*function_evaluations(:, tm);
    end
end
