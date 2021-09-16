function [xhat,outstat, XHAT] = runstreamestimate(dataout, m, options, Tf, k, eta, gamma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



%% Parse data and set all variables
XHAT={};
timevec = dataout.timevec;

X0 = ones(m,1); % initial guess
xhat = [];

%initialize variables
funarr = {};
% grdarr = {};
hesarr = {};

Kt_tvec =[] ;% k_t( ) time points 
Ktm1_tvec = [];
Kt_indvec = [] ;% k_t( ) indices 
% Ktm1_indvec = [] ;

Tdatind_t = [];% t-th Data frame indices 
% Tdatind_tm1 = []; 
Tdat_t = [];% t-th Data frame time points
% Tdat_tm1 = []; 
outstat=[];
for t=1:1:Tf
    % We curently maintain latency of at least frmlen/2 samples
     
    % Update frames 
    if t==1 %first block is half size and depends only on x1 f_1(x_1)
        % update data frames indices 
        Tdatind_t = (1:m/2);    
        Tdat_t = timevec(Tdatind_t);
        % Update basis functions (kernel) indices 
        Kt_indvec = 1:m;
        Kt_tvec = timevec(Kt_indvec);
    else
        % update data frames indices
        Tdatind_tm1 = Tdatind_t;
        Tdatind_t = Tdatind_tm1(end) + (1:m);
%         Tdat_tm1 = Tdat_t;
        Tdat_t = timevec(Tdatind_t);
        
        % Update basis functions (kernel) indices
        Ktm1_tvec = Kt_tvec;
        Ktm1_indvec = Kt_indvec;
        Kt_indvec = Ktm1_indvec(end) + (1:m);
        Kt_tvec =  timevec(Kt_indvec);
    end
            
    % Update bins spikes cound bt(i) = #spikes in [Tdat_t(i-1) , Tdat_t(i)]
    bt = dataout.spikevec(Tdatind_t);
    bt = bt(:);
    
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
    glbltvec = timevec(1:t*m);
    [reval, rgrad, rhes] = rkhsnorm(k, glbltvec, gamma);
%% Send to solver   
    % prep solver subroutine
    fun = @(X) Jcomb(X, J, Jgrad, Jhes, reval, rgrad, rhes, eta);
    
    % update guess
    Xc = [xhat ; X0];    
    sprintf('Runnig solver, t = % i\n',t)
    [xhat,fval,exitflag,output, grad, hessian ]  = fminunc(fun, Xc, options); %, [],[],[],[],Xc*0,[]);
    XHAT{t} = xhat;
    outstat(t,:) = norm(grad)/sqrt(t*m); %plot scaled error magnitude
        
end

end

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