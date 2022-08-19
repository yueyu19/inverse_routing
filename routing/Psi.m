function [ out ] = Psi(x, xref, str)
% compute the value and the gradient of the KL divergence between x and xref
% x = x + 1e-8; 
% xref = xref + 1e-8; 
% if strcmp(str, 'fval')
%     out = max(x'* (log(x)-log(xref)), 1e-16);
% elseif strcmp(str, 'grad')
%     out = log(x) - log(xref);
% end


if strcmp(str, 'fval')
    out = 0.5*sum((x-xref).^2);
elseif strcmp(str, 'grad')
    out = x-xref;
end

end

