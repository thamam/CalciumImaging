clear all;
clc
%%
m=6;
c=m/2+1;
NT = 4;
syms Kst [2*m ,2*m]
syms Kt   [m ,m]
syms delta
syms btm1 [m 1]
syms bt   [m 1]
syms xtm1 [m 1]
syms xt   [m 1]
syms jj
%tm1 => 1
% t => 2

Xt = [xtm1;xt];

temp = 0;
for i=c:(m+c-1)
    temp = temp + exp(delta*Kst(:,i).'*Xt);
end
ft = temp - bt.'*Kst(:,c:(m+c-1)).'*Xt;

grdf =gradient(ft,Xt);

hesf = hessian(ft, Xt);


%%
ytm1 = randn(m,1);
yt = randn(m,1);
tvec = linspace(0,2,2*m);
datavec = tvec(c:(m+c-1));
deltaval = tvec(2)-tvec(1);
btval = randi([1 3],m,1);
% Kt = kt(df_t(:)); % mx2m matrix
sig_f = 1;
l = 1;
k = @(t1,t2) sig_f * exp(-(t1-t2.').^2/(l^2));

K = k(tvec(:),tvec(:));
KK = K(:,c:m+c-1).';

onevec = ones(m,1);
fteval = @(xtm1, xt) onevec.'*exp(deltaval* KK*[xtm1;xt]) ...
    - btval.'*(KK*[xtm1;xt]);
ft_grad = @(xtm1,xt) KK.' *(deltaval * exp(deltaval* KK*[xtm1;xt])    - btval);
ft_hes = @(xtm1, xt) deltaval^2*KK.'*diag(exp(deltaval* KK*[xtm1;xt]))*KK;

%% verify ft - check
ftsymtemp = subs(ft,[delta ;xtm1; xt], [deltaval; ytm1;yt]) ;
ftsymtemp2 = subs(ftsymtemp,Kst,K);
ftsym = double(subs(ftsymtemp2, bt, btval)); 
ftnum = double(fteval(ytm1,yt));

diff_ft =  ftnum - ftsym;

%% verify gradient - check
ft_gradnum = double(ft_grad(ytm1,yt));
fgrdsymval = double(subs(grdf,[delta ;xtm1; xt; Kst(:); bt] , [deltaval; ytm1; yt; K(:) ;btval]));
grad_diff = ft_gradnum - fgrdsymval

%% Verify Hessian

ft_hesnum= double(ft_hes(ytm1,yt));
fhesSymval = double(subs(hesf ,[delta ;xtm1; xt; Kst(:); bt] , [deltaval; ytm1; yt; K(:) ;btval]));
hes_diff = ft_hesnum - fhesSymval






