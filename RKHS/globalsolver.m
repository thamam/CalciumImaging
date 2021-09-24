function [xhat,fval, exitflag] = globalsolver(X0, dataout, ktvec, k, eta)

N = numel(X0);
[f,g,H] = globalcostfunction(dataout, ktvec, k, N, eta);

alg2use = 'Quasi-Newton' ;
options = optimoptions('fminunc','Algorithm',alg2use,'SpecifyObjectiveGradient',true);
options.HessianFcn = 'objective';

fun = @(X) Jcomb(X,f,g,H);
[xhat, fval, exitflag ]  = fminunc(fun, X0, options); %, [],[],[],[],Xc*0,[]);

end


function [f, g, H] = globalcostfunction(dataout, ktvec, k, N, eta)
%GLOBALSOLVER Summary of this function goes here
%   Detailed explanation goes here
% N - number of basis function - points
tspikes = dataout.tspikes;
delta = dataout.delta;
%%%%
%
Mt = numel(tspikes); %number of events in data frame
% verify  sizes of kenerls returned values
% assert(size(ktm1(1),1)==m);
ut = @(t) k(ktvec(:), t);  % (2mx1)
wt = @(t) k(ktvec(:), t);  % (2mx1) %same as u_i but for tau(spikes)
Ut = ut(ktvec(:).');
Wt = wt(tspikes(:).');
v = @(X) exp(Ut.'*X);
%ft
f = @(X) ones(N,1).'*delta*v(X)- ones(Mt,1).'*Wt.'*X +eta*X.'*Ut*X ;
if nargout > 1 % gradient required
    % \nabla ft
    g = @(X) delta * Ut*v(X) - Wt*ones(Mt,1) + 2*eta*Ut*X;
    if nargout > 2 % Hessian required        
        % \nabla^2 ft
        Sigma = @(X) delta*diag(v(X));
        H =@(X) Ut * Sigma(X) * Ut.' + 2*eta*Ut;
    end
end
end

function [f, g, H] = Jcomb(X, f, g, H)
% Calculate objective f
f = f(X);
if nargout > 1 % gradient required
    g = g(X);
    if nargout > 2 % Hessian required
        H = H(X);
    end
end
end



