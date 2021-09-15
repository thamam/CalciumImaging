function [HesJ] = assemglobHes(hesarr, t, m, X)
%assemglobobj build up the global objective
%   Detailed explanation goes here
Imxm = eye(m*t);
Ui_ = @(ii) Imxm(:, (ii-2)*m + (1:2*m));
x_ = @(t) X((1:m) + (t-1)*m) ; 
U1 = Imxm(:,1:m);
HesJ =  U1*hesarr{1}(x_(1))*(U1.');
for ii=2:t      
    %% lift into apropriate space
    Ui = Ui_(ii);   
    HesJ =   HesJ +  Ui*hesarr{ii}(x_( ii-1), x_( ii))*(Ui.');        
end


