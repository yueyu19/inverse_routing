clear;clc;

pa = struct; % structure var for storing all parameters
pa.p = 3; % # of players
pa.E = [ 1, -1,  0,  0,  1, -1;
        -1,  1,  1, -1,  0,  0;
         0,  0, -1,  1, -1,  1]; % incidence matrix
  
[pa.n, pa.m] = size(pa.E); % # of nodes and links in the graph 
  

pa.b = 0.1*ones(pa.p*pa.m, 1);
pa.C = zeros(pa.p*pa.m, pa.p*pa.m);
pa.lam = 0.01;
pa.smat = [ 1,  0, -1;
           -1,  1,  0;
            0, -1   1];
pa.s = pa.smat(:);

opts = optimoptions('fsolve', 'MaxIterations', 1e3, 'MaxFunctionEvaluations', 1e8);
%%

[x, v, fval, output] = myfun(pa, zeros(pa.p*pa.m, 1), zeros(pa.p*pa.n, 1), opts);
%reshape(x, [], pa.p)

xref = [0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0]';

alp = 1e-2;
maxiter = 10;
Cplus = pa.C;
Psival = zeros(1, maxiter);
rho = 10;

for iter = 1:maxiter
    pa.C = Cplus;
    [x, v, fval, output] = myfun(pa, x, v, opts);
    max(abs(fval))
    if fval >1e-4
        break
    end
    gPsi = Psi(x, xref, 'grad');
    Cplus = projD (pa.C - alp*grad(x, v, gPsi, pa), rho);
    Psival(iter) = Psi(x, xref, 'fval');
end

figure
plot(1:maxiter, Psival, 'r', 'LineWidth', 2)






