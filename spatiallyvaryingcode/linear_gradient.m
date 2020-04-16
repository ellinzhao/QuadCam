function [f, g] = linear_gradient(x, A, Ah, b)
    u = A(x);
    r = u-b;
    g = Ah(r);
    %g = g/norm(g,'fro')*10000;
    f = sum(r(:).^2);
    %f = norm(r,'fro');