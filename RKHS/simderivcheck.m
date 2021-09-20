clear all;
clc
%
m=4; % #basis functions / frame
M = 3; %#of events in dataframe (actual spikes times, not restrcited to grid)
d = m; %length of dataframe (in kenel units
c=m/2+1;   NT = 4;
syms xtm1 [m 1]
syms xt   [m 1]
syms tau [ M 1]
syms tvec [d, 1]
delta = 0.5;
syms tvecktm1 [m, 1]
syms tveckt [m, 1]


% t=10;
t=1; %% uncomment to tst for the case ft = f1 (shorter block supported only on x1)


% define function handles

sig_f = 0.5;
sig_l = 1.5;
k = @(t1,t2) sig_f * exp(-(t1-t2).^2/(sig_l^2));

ai = @(t) k(tvecktm1, t); %(mx1)
bi = @(t) k(tveckt, t);   %(mx1)
X = [xtm1;xt];              %(2mx1)
% vi = @(t) [ai(t).' , bi(t).'];  % (1x2m)
ui = @(t) [ai(t) ; bi(t)];  % (2mx1)
ci = @(t) k(tvecktm1, t); %(mx1)
di = @(t) k(tveckt, t); %(mx1)
wi = @(t) [ci(t) ; di(t)];  % (2mx1) %same as u_i but for tau(spikes)
%% Let the computer compute the derivtives
ftsym=0;
for ii = 1:d
    ftsym = ftsym + delta*exp(ui(tvec(ii)).'*X)  ;
    %     ftsym = ftsym + delta*exp(vi(tvec(ii))*X)  ;
end

for mm = 1:M
    ftsym = ftsym - ci(tau(mm)).'*xtm1 - di(tau(mm)).'*xt;
end
% matrix form
ftsymMat = ones(m,1).'*delta*exp(ui(tvec.').'*X) - ones(M,1).'*wi(tau.').'*X;
% simplify(ftsym - ftsymMat)

grdf =gradient(ftsym,X);

hesf = hessian(ftsym, X);


%% Compute analytical derivatives
% ui = @(t) [ai(t) ; bi(t)];  % (2mx1)
% wi = @(t) [ci(t) ; di(t)];  % (2mx1) %same as u_i but for tau(spikes)
if t==1
%     X = subs(X,xtm1,zeros(m,1));
    x1=xt;
    U1 = bi(tvec.');
    W1 = di(tau.');
    v1 = exp(bi(tvec.').'*x1);    

    %ft
    numft = ones(m,1).'*delta*v1- ones(M,1).'*W1.'*x1;
    % \nabla ft
    numgrdf = delta * U1*v1 - W1*ones(M,1);
    
    % \nabla^2 ft
    Sigma1 = delta*diag(v1);
    numHesf =U1 * Sigma1 * U1.';
else
    Ui = ui(tvec.');
    Wi = wi(tau.');
    v = exp(ui(tvec.').'*X);
    %ft
    numft = ones(m,1).'*delta*v- ones(M,1).'*Wi.'*X;
    % \nabla ft
    numgrdf = delta * Ui*v - Wi*ones(M,1);
    
    % \nabla^2 ft
    Sigma = delta*diag(v);
    numHesf =Ui * Sigma * Ui.';
end
%% Compare symbolic and analytical derivative
% generate some random values to replace in expressions
if t==1
    ytm1 = zeros(m,1); 
else       
    ytm1 = randn(m,1);
end
yt = randn(m,1);
tauval  = randsample(linspace(0,2,100),M);
tauval=tauval(:);
tvec2m = linspace(0,2,2*m).';
tvecktm1val = tvec2m(1:m); %linspace(0,1,m).';
tvecktval =  tvec2m(m+1:2*m); %linspace(1,2,m).';
tvecval = tvec2m((1:m)+round(m/2));

%
%%%%%%%%%
%
ft_diff = (numft - ftsym);
subsft_diff = double(subs(ft_diff, [xtm1; xt; tau ; tvec; tvecktm1; tveckt], ...
    [ytm1; yt; tauval; tvecval ; tvecktm1val; tvecktval]));

if t==1    
    gradft_diff = (grdf(m+1:end) - numgrdf);
    % ft_diff_subed = subs(ft_diff,[xtm1; xt; tau ; tvec; tvecktm1; tveckt], ...
    %     [ytm1; yt; tauval; tvecval ; tvecktm1val; tvecktval])
    % double(ft_diff_subed)
    subsgrad = double(subs(gradft_diff, [xtm1; xt; tau ; tvec; tvecktm1; tveckt], ...
        [ytm1; yt; tauval; tvecval ; tvecktm1val; tvecktval]));
    
    hes_diff = (hesf(m+1:end,m +1:end) - numHesf );
    
    subsHesdiff = double(subs(hes_diff, [xtm1; xt; tau ; tvec; tvecktm1; tveckt], ...
        [ytm1; yt; tauval; tvecval ; tvecktm1val; tvecktval]));
else
    gradft_diff = (grdf - numgrdf);
    % ft_diff_subed = subs(ft_diff,[xtm1; xt; tau ; tvec; tvecktm1; tveckt], ...
    %     [ytm1; yt; tauval; tvecval ; tvecktm1val; tvecktval])
    % double(ft_diff_subed)
    subsgrad = double(subs(gradft_diff, [xtm1; xt; tau ; tvec; tvecktm1; tveckt], ...
        [ytm1; yt; tauval; tvecval ; tvecktm1val; tvecktval]));
    
    hes_diff = (numHesf - hesf);
    
    subsHesdiff = double(subs(hes_diff, [xtm1; xt; tau ; tvec; tvecktm1; tveckt], ...
        [ytm1; yt; tauval; tvecval ; tvecktm1val; tvecktval]));
end
%%

% Itm1 = [1 2];
% It   = [2,3];
% syms Kst [2*m ,2*m]
% syms Kt   [m ,m]
% syms kt [m 1]
% syms ktm1 [m 1]
% syms btm1 [m 1]
% syms bt   [m 1]
% % % % %
% % % %
% % % % datavec = tvec(c:(m+c-1));
% % % % deltaval = tvec(2)-tvec(1);
% % % % btval = randi([1 3],m,1);
% % % %
% % % %
% % % % K = k(tvec(:),tvec(:));
% % % % KK = K(:,c:m+c-1).';
% % % %
% % % % onevec = ones(m,1);
% % % % fteval = @(xtm1, xt) onevec.'*exp(deltaval* KK*[xtm1;xt]) ...
% % % %     - btval.'*(KK*[xtm1;xt]);
% % % % ft_grad = @(xtm1,xt) KK.' *(deltaval * exp(deltaval* KK*[xtm1;xt])    - btval);
% % % % ft_hes = @(xtm1, xt) deltaval^2*KK.'*diag(exp(deltaval* KK*[xtm1;xt]))*KK;
% % % %
% % % % % verify ft - check
% % % % ftsymtemp = subs(ft,[delta ;xtm1; xt], [deltaval; ytm1;yt]) ;
% % % % ftsymtemp2 = subs(ftsymtemp,Kst,K);
% % % % ftsym = double(subs(ftsymtemp2, bt, btval));
% % % % ftnum = double(fteval(ytm1,yt));
% % % %
% % % % diff_ft =  ftnum - ftsym;
% % % %
% % % % % verify gradient - check
% % % % ft_gradnum = double(ft_grad(ytm1,yt));
% % % % fgrdsymval = double(subs(grdf,[delta ;xtm1; xt; Kst(:); bt] , [deltaval; ytm1; yt; K(:) ;btval]));
% % % % grad_diff = ft_gradnum - fgrdsymval
% % % %
% % % % % Verify Hessian
% % % %
% % % % ft_hesnum= double(ft_hes(ytm1,yt));
% % % % fhesSymval = double(subs(hesf ,[delta ;xtm1; xt; Kst(:); bt] , [deltaval; ytm1; yt; K(:) ;btval]));
% % % % hes_diff = ft_hesnum - fhesSymval
% % % %
% % % %




