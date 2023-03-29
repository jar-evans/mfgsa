function U = norm_samples(Z)

    U = zeros(size(Z));
    input_ranges = [0.25 1.25 0.5];
    input_means = [0.5 25/8 -1];

    Z = Z - input_means;
    U(:,1) = Z(:,1) / input_ranges(1) /2;
    U(:,2) = Z(:,2) / input_ranges(2) /2;
    U(:,3) = Z(:,3) / input_ranges(3) /2;



end

