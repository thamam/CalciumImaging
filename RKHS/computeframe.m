function [unepsind, mintime] = computeframe(eps, delta, k )
%computeframe
% Find the minimal distance where k(t,t+delta_t)<eps;
%
%%

% eps = 1e-3;
% k = @(i,j) sig_f*exp(-(i-j).^2./sig_l^2);
% delta = 30e-3;
% sig_f = 1;
% sig_l= 1/2;;
% w = @(i,j) -(i-j).^2./sig_l^2;

N = 100;
for i=1:10
    timebrackets = (0:delta:N*delta);
    keval = k(0,timebrackets);
%     plot(timebrackets,keval);
    unepsind = find(keval<eps,1,'first');
    if ~isempty(unepsind)
        mintime = unepsind*delta;
        return;
    else
        N=N*2;
    end
end

 if isempty(unepsind)
     error('can not find k<eps in range\n');
 end


end

