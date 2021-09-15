function [reval, rgrad, rhes] = rkhsnorm(k, glbltvec, gamma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('gamma','var')
    gamma = 0;
end
    

K = k(glbltvec(:), glbltvec(:).') + gamma*eye(numel(glbltvec));

reval = @(X) X.' * K * X ;

rgrad =@(X) 2*K * X;

rhes = 2*K;



end

