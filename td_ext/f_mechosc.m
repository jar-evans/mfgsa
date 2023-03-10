function [y, t] = f_mechosc(I, T, dt)

    a = I(:, 1);
    b = I(:, 2);
    l = I(:, 3);
    t = 0:dt:T;

    y = l.*exp(-a*t).*(cos(b*t) + (a/b)*sin(b*t));
end