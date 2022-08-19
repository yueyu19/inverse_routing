function [x, v, fval, output] = myfun(pa, x0, v0, opts)

% x is the Nash defined by pa

[equ, fval, ~, output] = fsolve(@nash, [x0; v0], opts);

x = equ(1:pa.p*pa.m);
v = equ(pa.p*pa.m+1:end);

function y = nash(xi)

y = [xi(1:pa.p*pa.m) - exp((1/pa.lam)*(kron(eye(pa.p), pa.E')*xi(pa.p*pa.m+1:end)-pa.b-pa.C*xi(1:pa.p*pa.m))-ones(pa.p*pa.m, 1));
     pa.s - kron(eye(pa.p), pa.E)*xi(1:pa.p*pa.m)];
end

end