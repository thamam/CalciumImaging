function [f, g, H] = computelossatevetns(delta, dttvec, tau_t, ktm1, kt, m, t )
%COMPUTELOSSATEVETNS Summary of this function goes here
%   Loss function using interarrival approach. i.e., in the second ML term we use s(\tau_m)
% Inputs: (delta, dttvec, tau_t, ktm1,kt,m,t)
% % X = [xtm1;xt];

%%%%
%
Mt = numel(tau_t); %number of events in data frame
d = numel(dttvec);
% verify  sizes of kenerls returned values
% assert(size(ktm1(1),1)==m); 
assert(size(kt(1),1)==m);
dttvec = dttvec(:);
tau_t = tau_t(:);

at = ktm1;  %(mx1)
bt = kt;   %(mx1)
ut = @(t) [at(t) ; bt(t)];  % (2mx1)

ct = ktm1; %@(t) k(tvecktm1, t);
dt = kt; % @(t) k(tveckt, t);
wt = @(t) [ct(t) ; dt(t)];  % (2mx1) %same as u_i but for tau(spikes)


if t==1 %shoter block f1(x) vs ft(xtm1,xt)
    U1 = bt(dttvec.');
    W1 = dt(tau_t.');
    v1 = @(X) exp(bt(dttvec.').'*X);
    
    %ft
    f = @(x1) ones(d,1).'*delta*v1(x1)- ones(Mt,1).'*W1.'*x1;
    % \nabla ft
    g = @(x1) delta * U1*v1(x1) - W1*ones(Mt,1);
    
    % \nabla^2 ft
    Sigma1 = @(X) delta*diag(v1(X));
    H = @(x1) U1 * Sigma1(x1) * U1.';
else
    Ut = ut(dttvec.');
    Wt = wt(tau_t.');
    v = @(X) exp(ut(dttvec.').'*X);
    
    %ft
    f = @(xtm1, xt) ones(d,1).'*delta*v([xtm1;xt])- ones(Mt,1).'*Wt.'*[xtm1;xt];
    % \nabla ft
    g = @(xtm1, xt) delta * Ut*v([xtm1;xt]) - Wt*ones(Mt,1);
    
    % \nabla^2 ft
    Sigma = @(X) delta*diag(v(X));
    H =@(xtm1, xt) Ut * Sigma([xtm1;xt]) * Ut.';
end
end

