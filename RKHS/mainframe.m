close all;
clear;
clc;
eps = 1e-3;
% delta = 30e-3;
sig_f = 1;
sig_l= 1/2;
k = @(i,j) sig_f*exp(-(i-j).^2./sig_l^2);

eta = 2; % penalty term wieght
gamma = 1e-6; %Tikhonov regularization to stabilize Kernel matrix
options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
options.HessianFcn = 'objective';
% options.CheckGradients = true;

options.FiniteDifferenceType = 'central';
options.FiniteDifferenceStepSize = 1e-10;
options.OptimalityTolerance = 1e-10;
options.MaxIterations = 1000;
options.FunctionTolerance = 1e-10;
%% read data

dataname = 'test';
[dataout] = loaddata(dataname) ; %read data for code testing 

[cind, mintime] = computeframe(eps, dataout.delta, k );

frmlen = 2*mintime;
m = 2*cind; %frame length / # of basis functions

Tf = floor(dataout.tf/frmlen); %adjust sim length to full length frames
X0 = ones(m,1);
basisframeindMat = [];
basisframetimeMat = [];
funarr = {};
grdarr = {};
hesarr = {};
%% Streaming loop
Xc=X0;
for t=1:1:Tf
    % We curently maintain latency of at least 45 samples (~45*delta[s])
 
    %% compute ft and derivatives 
    if t==1 %first block is half size
        timeind_t = (t-1)*(m/2)+(1:m/2);        
    else
        timeind_t = m/2 + (t-2)*(m)+(1:m);
    end
    %xt(baasis time array)
    basisframeind = (1:m) + (t-1)*m;
    basisframetime = dataout.deltatiks(basisframeind);
    basisframeindMat(:,t) = basisframeind(:);
    basisframetimeMat(:,t) = basisframetime(:);
    %compute data frame tiks 
    df_t = dataout.deltatiks(timeind_t);
    %Extract spikes at current frame
    bt = dataout.spikevec(timeind_t);
    %compute basis functions (corresponds to Xt) indices and update    
    if t==1
        kt = @(tj) k(basisframetime(:).', tj);
        [fteval,ft_grad, ft_hes] = ...
        buildft(dataout.delta, df_t, bt, [], kt, m, t);        
    else        
        ktm1 = @(tj) k(basisframetimeMat(:,t-1).', tj);
        kt = @(tj) k(basisframetime(:).', tj);                
        [fteval,ft_grad, ft_hes] = ...
        buildft(dataout.delta, df_t, bt, ktm1,kt,m,t);
    end    
    funarr{t} = fteval ;
    garr{t} = ft_grad ;  
    hesarr{t} = ft_hes ;
    %% update objective and derivatives
    if t==1
        J = fteval;
        gJ = ft_grad;
        HesJ = ft_hes;
    else
        [J] = @(X) assemglobobj(funarr, t, m, X);
        [gJ] = @(X) assemglobgrad(garr, t, m, X);
        [HesJ] =@(X) assemglobHes(hesarr, t, m, X);    
    end
    % prepare rkhs reg term
    glbltvec = dataout.deltatiks(1:t*m);
    [reval, rgrad, rhes] = rkhsnorm(k, glbltvec, gamma);
    %% Send to solver   
    fun = @(X) Jcomb(X, J, gJ, HesJ, reval, rgrad, rhes, eta);
    [xhat,fval,exitflag,output, grad, hessian ]  = fminunc(fun, Xc, options); %, [],[],[],[],Xc*0,[]);
    norm(grad)/sqrt(t*m)
    Xc = [xhat ; X0];
        
end

%% Plot result
% Compute lambdaVec
tplot = dataout.ts:(dataout.delta/100):dataout.tf;
svec = k(tplot(:),glbltvec(:).')*xhat;
lamvec = exp(svec);
%%
figure(1),clf
plot(tplot,lamvec,'--b')
hold all
stem(dataout.deltatiks(1:t*m),5*dataout.spikevec((1:t*m)))



function [f, g, H] = Jcomb(X, J, gJ, HesJ, reval, rgrad, rhes, eta)
% Calculate objective f
f = J(X) + eta*reval(X);

if nargout > 1 % gradient required
    g = gJ(X)+eta*rgrad(X);
    
    if nargout > 2 % Hessian required
        
        H = HesJ(X)+eta*rhes;
    end
end

end