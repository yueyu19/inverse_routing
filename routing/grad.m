function [gC] = grad(x, v, gPsi, pa)

D = diag(exp((1/pa.lam)*(kron(eye(pa.p), pa.E')*v-pa.b-pa.C*x)-ones(pa.p*pa.m, 1)));
J = [eye(pa.p*pa.m)+(1/pa.lam)*D*pa.C, (1/pa.lam)*D*kron(eye(pa.p), pa.E');
     -kron(eye(pa.p), pa.E), zeros(pa.p*pa.n, pa.p*pa.n)];
 
gC = -(1/pa.lam)*[D', zeros(pa.p*pa.m, pa.p*pa.n)]*pinv(J)'*[gPsi; zeros(pa.p*pa.n, 1)]*x'; 

end

