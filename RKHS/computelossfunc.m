function [fteval,ft_grad, ft_hes] = computelossfunc(delta, Tdat_t, bt, ktm1,kt,m,t)
%buildft prepare ft function handle
%   Detailed explanation goes here


if t==1 %shoter block f1(x) vs ft(xtm1,xt)
   
    onevec = ones(m/2,1);

    Kt = kt(Tdat_t(:)); % mx2m matrix
    
    fteval = @(xt) delta*onevec.'*exp( Kt*xt) ...
        - bt.'*(Kt*xt);
        
    ft_grad = @(xt) delta * Kt.' * exp( Kt*xt)    - bt);
    
    
    ft_hes = @(xt) delta *(Kt.'*exp(Kt*xt)*exp( Kt*xt).'*Kt);
    
else    % t=2,3,...
    onevec = ones(m,1);
    
    Kt = [ktm1(Tdat_t(:)) , kt(Tdat_t(:)) ]; % mx2m matrix    
    
    fteval = @(xtm1, xt) delta*onevec.'*exp( Kt*[xtm1;xt]) ...
        - bt.'*(Kt*[xtm1;xt]);  
    
    ft_grad = @(xtm1,xt) delta*Kt.' *( exp( Kt*[xtm1;xt])    - bt);
        
    ft_hes = @(xtm1, xt) delta*(Kt.'*diag(exp( Kt*[xtm1;xt]))*Kt);    
end


end

