function [xhat,outstat,ktvec, dtvec, XHAT, tauArray] = ...
    runstreamestimate(dataout, m, options, Tf, k, eta, gamma, BUFFERLENGTH)

% runstreamestimate Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse data and set all variables

XHAT    = zeros(Tf*m,Tf);                                                  % 
timevec = dataout.timevec;                                                 % 
tspikes = dataout.tspikes;                                                 % 
X0      = ones(m,1);                                                       % initial guess
xhat    = [];                                                              % 
% initialize variables
funarr   = {};                                                             % 
hesarr   = {};                                                             % 
kttvec   = [];                                                             % k_t( ) time points
ktm1tvec = [];                                                             % 
ktind    = [];                                                             % k_t( ) indices
% grdarr = {};
% Ktm1_indvec = [] ;

dtind    = [];                                                             % t-th Data frame indices
dttvec   = [];                                                             % t-th Data frame time points
outstat  = [];                                                             % 
tauArray = {};                                                             % 
% Tdatind_tm1 = [];
% Tdat_tm1 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start stremaing

for t = 1:1:Tf                                                             % Loop over batches 
    % We curently maintain latency of at least frmlen/2 samples
    % Update frames
    if t==1                                                                % The first block is half size and depends only on x1 f_1(x_1)
        % update data frames indices
        dtind  = (1:m/2);
        dttvec = timevec(dtind);
        % Update basis functions (kernel) indices
        ktind  = 1:m;
        kttvec = timevec(ktind);
    else
        % update data frames indices
        Tdatind_tm1 = dtind;
        dtind = Tdatind_tm1(end) + (1:m);
        %         Tdat_tm1 = Tdat_t;
        dttvec = timevec(dtind);   
        % Update basis functions (kernel) indices
        ktm1tvec = kttvec;
        Ktm1_indvec = ktind;
        ktind = Ktm1_indvec(end) + (1:m);
        kttvec =  timevec(ktind);
    end
    % Update bins spikes cound bt(i) = #spikes in [Tdat_t(i-1) , Tdat_t(i)]
    tframe      = [dttvec(1),dttvec(end)];
    tau_t       = tspikes(tspikes>=tframe(1) & tspikes<=tframe(2));
    tauArray{t} = tau_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Update basis functions frame handles and compute ft handle
    if t==1
        kt = @(tj) k(kttvec(:), tj);
        [fteval,ft_grad, ft_hes] = computelossatevetns(dataout.delta, dttvec, tau_t, [], kt, m, t);
    else
        ktm1 = @(tj) k(ktm1tvec(:), tj);
        kt   = @(tj) k(kttvec(:), tj);
        [fteval,ft_grad, ft_hes] = computelossatevetns(dataout.delta, dttvec, tau_t, ktm1, kt, m, t);
    end
    % Add ft and derivatives to storage array
    funarr{t} = fteval ;
    garr{t}   = ft_grad ;
    hesarr{t} = ft_hes ;    
    %% update global objective and derivatives
    istart = max(t-BUFFERLENGTH+1,1);  
    if t==1
        J     = fteval;
        Jgrad = ft_grad;
        Jhes  = ft_hes;
    else
        if istart==1 %buffer isn't active
            [J]     = @(X) assemglobobj(funarr, t, m, X );
            [Jgrad] = @(X) assemglobgrad(garr, t, m, X );
            [Jhes]  = @(X) assemglobHes(hesarr, t, m, X );
        else % assemble with truncated arrays
            [J]     = @(X) assemglobobj({{},funarr{istart:t}}, t, m, X, BUFFERLENGTH);
            [Jgrad] = @(X) assemglobgrad({{},garr{istart:t}}, t, m, X, BUFFERLENGTH );
            [Jhes]  = @(X) assemglobHes({{},hesarr{istart:t}}, t, m, X, BUFFERLENGTH );
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% prepare rkhs reg term
    % The reg term is treated as reg to the global problem instead of part
    % of the local loss to simplify approach
    glbltvec = timevec(1:t*m);
    [reval, rgrad, rhes] = rkhsnorm(k, glbltvec, gamma, BUFFERLENGTH, t, m, istart);
    %% Send to solver
    % prep solver subroutine
    
    if (istart>=2)
        Xtmb   = xhat(1:m);                                                % update x_{T-b} - latest inactive variable
        Xc     = [xhat(m+1:end) ; X0];                                     % update guess        
        funbuf = @(X) Jcomb([Xtmb;X],J,Jgrad,Jhes,reval,rgrad,rhes,eta);   % Prep opt call subroutine
        sprintf('Runnig solver, t = % i\n',t)
        [xhat, fval, exitflag, output, grad, hessian] = ...
                                             fminunc(funbuf, Xc, options); % Use fmincon to minimize the local cost  , [],[],[],[],Xc*0,[]);
        XHAT( (t-BUFFERLENGTH)*m+1:t*m,t) = xhat;                          % Save the estimate
    else
        Xc  = [xhat ; X0];
        fun = @(X) Jcomb(X, J, Jgrad, Jhes, reval, rgrad, rhes, eta);
        sprintf('Runnig solver, t = % i\n',t)
        [xhat, fval, exitflag, output, grad, hessian ] = ...
                                                fminunc(fun, Xc, options); %, [],[],[],[],Xc*0,[]);
        XHAT(1:t*m,t) = xhat;
    end
    
    
    outstat(t,:) = norm(grad)/sqrt(t*m); %plot scaled error magnitude
end
ktvec = timevec(1:ktind(end));
dtvec = timevec(1:dtind(end));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% Code snipets - delete once code is stable
% function [f, g, H] = JcombWbuff(Xtmb, X, J, gJ, HesJ, reval, rgrad, rhes, eta)
% % Calculate objective f
% f = J([Xtmb;X]) + eta*reval([Xtmb;X]);
% if nargout > 1 % gradient required
%     g = gJ([Xtmb;X])+eta*rgrad([Xtmb;X]);
%     if nargout > 2 % Hessian required
%         H = HesJ([Xtmb;X])+eta*rhes;
%     end
% end
% end

