function [ C ] = projD( C, rho )
% compute the projection of C onto the set of matrices, whose symmetric
% part is PSD

C1 = 0.5*(C-C');
C2 = 0.5*(C+C');

[V, D] = eig(C2);

C = real(C1 + V*max(D, 0)*V');

C = (rho/max(rho, norm(C, 'fro')))*C;

end

