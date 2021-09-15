close all;
clear;
clc;
%%
sig_f = 1;
l_w = .5;
n=10;

% syms a [n 1]
% syms b [n 1]

a = randn(n,1);
b = randi([0 1],n,1);
b=ones(n,1);

tbar = linspace(0,1,n+1).';
delta = tbar(2)-tbar(1);
tlamb = tbar(2:end) - delta/2;
%

k_d = @(d) sig_f*exp(-(d).^2./l_w^2);
v = @(t) k(tlamb, t);  % p(t)=psi(t) = [k(t_1,t), ..., k(t_n,t)]^T

t_mat = abs(tlamb - tlamb.');
K = k_d(t_mat);
rank(K)

%% computing P and it's first and second derivative in sum form 

%Penalty function
P_for = 0;
for i=1:n
     temp = delta*exp(a.'*v(tlamb(i))) - b(i)*a.'*v(tlamb(i)) ;
    P_for = P_for + temp;      
end
P_f = P_for + a.'*K*a;
% grad P_f
gradP_f=zeros(n,1);
for i=1:n
    temp = (delta*exp(a.'* v(tlamb(i)) )- b(i))*v(tlamb(i));
    gradP_f = gradP_f+temp;
end
gradP_f = gradP_f+2*K*a;
% Hess P_f
HesP_f=zeros(n,n);
for i=1:n
    temp = delta*exp(a.'* v(tlamb(i)))*v(tlamb(i))*v(tlamb(i)).';
    HesP_f = HesP_f+temp;
end
HesP_f = HesP_f+2*K;

%% computing P and it's first and second derivative in matrix form
r = exp(K*a);
del_1 = ones(n,1)*delta;
P_v = del_1.'*r -(K*a).'*b + a.'*K*a;

e_P = (P_v - P_f)/P_f;


%Grad P
z_op = delta*exp(K.'*a) - b;
gradP_v = K*z_op+2*K*a;

e_grad = norm(gradP_f-gradP_v)/norm(gradP_f);

% Hes p
dG = delta*diag(exp(K*a));
HesP_v=K * dG * K +2*K;

e_grad = norm(HesP_f-HesP_v)/norm(HesP_f);
%%
% syms A [3 3]
% A = randi(3,3);
% A = A+A.';
% syms x [3 1]
% C=K*x;
% f = exp(A*x);
% fgrad = jacobian(f)



%%
