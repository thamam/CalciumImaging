function [reval, rgrad, rhes] = rkhsnorm(k, glbltvec, gamma, BUFFERLENGTH, t, m, istart)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('gamma','var')
    gamma = 0;
end


if t-BUFFERLENGTH >0 %buffer is active
    glbltvec = glbltvec((istart-2)*m+1:end);
    K = k(glbltvec(:), glbltvec(:).') + gamma*eye(numel(glbltvec));
    
    reval = @(X) X.' * K * X;
   
    rgrad =@(X) 2*K(m+1:end,:) * X;
    rhes = 2*K(m+1:end,m+1:end);
else
    K = k(glbltvec(:), glbltvec(:).') + gamma*eye(numel(glbltvec));
    reval = @(X) X.' * K * X ;
    rgrad =@(X) 2*K * X;
    rhes = 2*K;
end


end

