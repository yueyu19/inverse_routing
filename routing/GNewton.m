function [ varargout ] = GNewton( x, pa )

nOutputs = nargout;
varargout = cell(1,nOutputs);

% Solve for the Nash equilibrium using the Newton method 

maxiter = 10000;
res = nan(maxiter, 1);
eps = 1e-4;

for iter = 1:maxiter
    
    
    [r, J] = eqnNash(x, pa);
    
    res(iter) = norm(r, Inf);
    
    if res(iter) < eps
        res(iter+1:end) = [];
        break
    end
    
    %x
    [d, ~] = lsqr(J, -r, 1e-4, 1e4); % compute Gauss-Newton descent direction
    
    % line search with Wolfe conditions
    
    c1 = 1e-4;
    c2 = 0.9;
    alp = 0;
    t = 1;
    bet = 1e4;
    
    
    for liter = 1:1000
        [rt, Jt] = eqnNash(x+t*d, pa);
%         fprintf('Newton iter %3d, line search iter %3d, step size %3.2e, error %3.2e \n', ...
%             iter, liter, t, norm(r))
       
        
        if 0.5*sum(rt.^2) > 0.5*sum(r.^2) + c1*t*r'*J*d
            bet = t;
            t = 0.5*(alp+bet);
        elseif rt'*Jt*d < c2*t*r'*J*d
            alp = t;
            if bet == 1e4
                t = 2*alp;
            else
                t = 0.5*(alp+bet);
            end
        else
            break
        end

        
        if max(t-alp, bet-t) < 1e-3
            break
        end
    end
    
    x = x + t*d;
end

% length(res)
% 
if length(res) == maxiter
    error('Maximum number of Newton iterations reached')
end

varargout{1} = x;

if nOutputs > 1
   varargout{2} = res;
end
    
end

