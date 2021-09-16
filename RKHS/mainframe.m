%% Streaming spikes regression - main frame
%mainframe.m
% written by Tomer H.Hamam tomerhamam@gmail.com 09.16.2021
% streaming estimation of intensity function from spikes measurements
% User options:
% Select data to process by setting the following parameters:
% dataname = 'test' / 'user' / 
% datapath
%
%%%%%
close all;
clear;
clc;
%% setting simulation and algorithm parameters
% Algorithm parameters
supeps = 1e-3; % support suppers threshold
sig_f = 1; sig_l= 1/2; %kernel initial parameters
k = @(i,j) sig_f*exp(-(i-j).^2./sig_l^2); %kernel handle function
eta = 2; % reg. penalty wieght
gamma = 1e-9; % Tikhonov regularization constant 

% Opt. solver params (for Matlab builtin solver)
alg2use = 'quasi-newton';
options = optimoptions('fminunc','Algorithm',alg2use,'SpecifyObjectiveGradient',true); 

% options.HessianFcn = 'objective';
% options.CheckGradients = true;
options.FiniteDifferenceType = 'central';
options.FiniteDifferenceStepSize = 1e-10;
options.OptimalityTolerance = 1e-10;
options.MaxIterations = 500;
options.FunctionTolerance = 1e-10;

%% Load data

% change data name into on
dataname = 'test'; 
[dataout] = loaddata(dataname) ; %read data for code testing 

% Compute frame size as 2*tmin, tmin = min_dt s.t. (k(|dt|)<supeps) use
% current sig_l, sig_f, delta check later if need to change frames' size
[frmlen, cind, mintime , m] = computeframe(supeps, dataout.delta, k ); 

% cind  - relative index to begining of basisframe  where dataframes
% starts\ends
% m - #basis functions/frame


% Set simulation length as N*frmlen N integer - truncate if needed 
Tf = floor(dataout.tf/frmlen); % simulation length
%% Streaming loop

X0 = ones(m,1); % initial guess
xhat = [];

%initialize variables
basisframeindMat = [];
basisframetimeMat = [];
funarr = {};
grdarr = {};
hesarr = {};

Kt_tvec =[] ;% k_t( ) time points 
Ktm1_tvec = [];
Kt_indvec = [] ;% k_t( ) indices 
Ktm1_indvec = [] ;

Tdatind_t = [];% t-th Data frame indices 
Tdatind_tm1 = []; 
Tdat_t = [];% t-th Data frame time points
Tdat_tm1 = []; 
for t=1:1:Tf
    % We curently maintain latency of at least frmlen/2 samples
     
    % Update frames 
    if t==1 %first block is half size and depends only on x1 f_1(x_1)
        % update data frames indices 
        Tdatind_t = (1:m/2);    
        Tdat_t = dataout.deltatiks(Tdatind_t);
        % Update basis functions (kernel) indices 
        Kt_indvec = 1:m;
        Kt_tvec =  dataout.deltatiks(Kt_indvec);
    else
        % update data frames indices
        Tdatind_tm1 = Tdatind_t;
        Tdatind_t = Tdatind_tm1(end) + (1:m);
        Tdat_tm1 = Tdat_t;
        Tdat_t = dataout.deltatiks(Tdatind_t);
        
        % Update basis functions (kernel) indices
        Ktm1_tvec = Kt_tvec;
        Ktm1_indvec = Kt_indvec;
        Kt_indvec = Ktm1_indvec(end) + (1:m);
        Kt_tvec =  dataout.deltatiks(Kt_indvec);
    end
            
    % Update bins spikes cound bt(i) = #spikes in [Tdat_t(i-1) , Tdat_t(i)]
    bt = dataout.spikevec(Tdatind_t);
    
    % Update basis functions frame handles and compute ft handle
    if t==1
        kt = @(tj) k(Kt_tvec(:).', tj);
        [fteval,ft_grad, ft_hes] = computelossfunc(dataout.delta, Tdat_t, bt, [], kt, m, t);        
    else        
        ktm1 = @(tj) k(Ktm1_tvec(:).', tj);
        kt   = @(tj) k(Kt_tvec(:).', tj);                
        [fteval,ft_grad, ft_hes] = computelossfunc(dataout.delta, Tdat_t, bt, ktm1,kt,m,t);
    end    
    
    % Add ft and derivatives to storage array
    funarr{t} = fteval ;
    garr{t} = ft_grad ;  
    hesarr{t} = ft_hes ;
    
%% update global objective and derivatives
    if t==1
        J = fteval;
        Jgrad = ft_grad;
        Jhes = ft_hes;
    else
        [J] = @(X) assemglobobj(funarr, t, m, X);
        [Jgrad] = @(X) assemglobgrad(garr, t, m, X);
        [Jhes] =@(X) assemglobHes(hesarr, t, m, X);    
    end
    
%% prepare rkhs reg term
    % The reg term is treated as reg to the global problem instead of part
    % of the local loss to simplify approach
    glbltvec = dataout.deltatiks(1:t*m);
    [reval, rgrad, rhes] = rkhsnorm(k, glbltvec, gamma);
%% Send to solver   

    % prep solver subroutine
    fun = @(X) Jcomb(X, J, Jgrad, Jhes, reval, rgrad, rhes, eta);
    
    % update guess
    Xc = [xhat ; X0];
    
    sprintf('Runnig solver, t = % i\n',t)
    [xhat,fval,exitflag,output, grad, hessian ]  = fminunc(fun, Xc, options); %, [],[],[],[],Xc*0,[]);
    
    norm(grad)/sqrt(t*m) %plot scaled error magnitude
        
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


%% subroutine for Matlab solver that returns both fun or its derivatives
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