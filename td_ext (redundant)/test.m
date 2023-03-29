
for III = 1:100
    III
    mechanical_oscillator;
    A = G_alg2(1, :);
    B = G_alg2(2, :);
    L = G_alg2(3, :);

    if (III == 1)
        a = A;
        b = B;
        l = L;
    else
        a = vertcat(a, A);
        b = vertcat(b, B);
        l = vertcat(l, L);
    end
end

plot(max(a)); hold on;
plot(min(a));
plot(max(b));
plot(min(b));
plot(max(l));
plot(min(l)); hold off;