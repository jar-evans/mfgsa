function y = f_mechosc(I, t)

    a = I(:, 1);
    b = I(:, 2);
    l = I(:, 3);

    y = l.*exp(-a*t).*(cos(b*t) + (a/b)*sin(b*t));
end