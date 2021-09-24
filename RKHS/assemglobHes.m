function [HesJ] = assemglobHes(hesarr, t, m, X, BUFFERSIZE)
%assemglobobj build up the global objective
%   When buffer is active, for 1st block- need to keep only part-deriv
%   w.r.t. x_{T-b+1}
% we do that by computing with it at HesJ_, and then discarding info
% relevant to x_{T-b}.Detailed explanation goes here
x_ = @(t) X((1:m) + (t-1)*m) ;
if exist('BUFFERSIZE','var') && (t-BUFFERSIZE)>0 %BUFFER is active
    Imxm = eye(m*(BUFFERSIZE+1));
    Ui_ = @(ii) Imxm(:, (ii-2)*m + (1:2*m));          
        
    HesJ_ =  zeros(m*(BUFFERSIZE+1), m*(BUFFERSIZE+1));
    
    for ii=2:BUFFERSIZE+1
        %% lift into apropriate space
        Ui = Ui_(ii);
        HesJ_ =   HesJ_ +  Ui*hesarr{ii}(x_( ii-1), x_( ii))*(Ui.');
    end
    HesJ = HesJ_(m+1:end,m+1:end);
    
else
     
    Imxm = eye(m*t);
    Ui_ = @(ii) Imxm(:, (ii-2)*m + (1:2*m));
    U1 = Imxm(:,1:m);
    HesJ =  U1*hesarr{1}(x_(1))*(U1.');
    for ii=2:t
        %% lift into apropriate space
        Ui = Ui_(ii);
        HesJ =   HesJ +  Ui*hesarr{ii}(x_( ii-1), x_( ii))*(Ui.');
    end
    
    
end

