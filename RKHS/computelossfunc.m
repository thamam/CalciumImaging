function [fteval,ft_grad, ft_hes] = computelossfunc(delta, Tdat_t, bt, ktm1,kt,m,t)
%buildft prepare ft function handle
%   Detailed explanation goes here


if t==1 %shoter block f1(x) vs ft(xtm1,xt)
   
    onevec = ones(m/2,1);

    Kt = kt(Tdat_t(:)); % mx2m matrix
    
    fteval = @(xt) onevec.'*exp(delta* Kt*xt) ...
        - bt.'*(Kt*xt);
        
    ft_grad = @(xt) Kt.' *(delta * exp(delta* Kt*xt)    - bt);
    
    
    ft_hes = @(xt) delta^2 *(Kt.'*exp(delta* Kt*xt)*exp(delta* Kt*xt).'*Kt);
    
else    % t=2,3,...
    onevec = ones(m,1);
    
    Kt = [ktm1(Tdat_t(:)) , kt(Tdat_t(:)) ]; % mx2m matrix    
    
    fteval = @(xtm1, xt) onevec.'*exp(delta* Kt*[xtm1;xt]) ...
        - bt.'*(Kt*[xtm1;xt]);  
    
    ft_grad = @(xtm1,xt) Kt.' *(delta * exp(delta* Kt*[xtm1;xt])    - bt);
        
    ft_hes = @(xtm1, xt) delta^2 *(Kt.'*diag(exp(delta* Kt*[xtm1;xt]))*Kt);    
end


end

