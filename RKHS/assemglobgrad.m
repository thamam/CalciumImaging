function [gJ] = assemglobgrad(garr, t, m, X, BUFFERSIZE)
%assemglobobj build up the global objective
    %for 1st block- need to keep only part-deriv  w.r.t. x_{T-b+1}
    % we do that by computing with it at gJ_, and then discarding first
    % block.

x_ = @(t) X((1:m) + (t-1)*m) ;

if exist('BUFFERSIZE','var') && (t-BUFFERSIZE)>0 %BUFFER is active
    Imxm = eye(m*(BUFFERSIZE+1));
    Ui_ = @(ii) Imxm(:, (ii-2)*m + (1:2*m));

    gJ_=zeros(m*(BUFFERSIZE+1),1);
    for ii=2:BUFFERSIZE+1
        Ui = Ui_(ii);
        gJ_ =  gJ_ + Ui*garr{ii}(x_( ii-1),x_( ii));
    end
    gJ = gJ_(m+1:end);
    
else   
    Imxm = eye(m*t);
    Ui_ = @(ii) Imxm(:, (ii-2)*m + (1:2*m));
    U1 = Imxm(:,1:m);
    gJ = U1*garr{1}(x_( 1));
    for ii=2:t %% lift into apropriate space
        Ui = Ui_(ii);
        gJ =  gJ + Ui*garr{ii}(x_( ii-1),x_( ii));
    end
end