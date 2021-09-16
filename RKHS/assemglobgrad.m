function [gJ] = assemglobgrad(garr, t, m, X)
%assemglobobj build up the global objective
%   Detailed explanation goes here
Imxm = eye(m*t);

Ui_ = @(ii) Imxm(:, (ii-2)*m + (1:2*m));

x_ = @(t) X((1:m) + (t-1)*m) ; 

U1 = Imxm(:,1:m);
gJ = U1*garr{1}(x_( 1));
for ii=2:t   
 
    
    %% lift into apropriate space
    Ui = Ui_(ii);

    gJ =  gJ + Ui*garr{ii}(x_( ii-1),x_( ii));
    
        
end
