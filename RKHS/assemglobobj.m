function [J] = assemglobobj(funarr, t, m, X, BUFFERSIZE)
%assemglobobj build up the global objective
%   Detailed explanation goes here

x_ = @(t) X((1:m) + (t-1)*m) ;

if exist('BUFFERSIZE','var') && (t-BUFFERSIZE)>0 %BUFFER is active
    J = 0;
%     funbuf = {{},funarr{:}}; %push by one to align with variables numbering
    for ii=2:BUFFERSIZE+1   
        J =  J + funarr{ii}(x_( ii-1),x_( ii));
    end
else
    J = funarr{1}(x_(1));
    
    for ii=2:t
        J =  J + funarr{ii}(x_( ii-1),x_( ii));
    end
end


% % % function [J, gJ, HesJ] = assemglobobj(funarr, garr , hesarr, t, m)
% % % assemglobobj build up the global objective
% % %   Detailed explanation goes here
% % % Imxm = eye(m*t);
% % %
% % % Ui_ = @(ii) Imxm((ii-1)*t + 1:m,:);
% % %
% % % x_ = @(X,t) X((1:m) + (t-1)*m) ;
% % %
% % % J = @(X) funarr{1}(x_(X, 1));
% % % U1 = Ui_(1);
% % % gJ = @(X) U1(1)*garr{1}(x_(X, 1));
% % % HesJ =  @(X) U1(1)*hesarr{1}(x_(X, 1))*(U1.');
% % % for i=2:t
% % %     J = @(X) J(X) + funarr{i}(x_(X, i-1),x_(X, i));
% % %
% % %     % lift into apropriate space
% % %     Ui = Ui_(i);
% % %
% % %     gJ = @(X) gJ(X) + Ui*garr{1}(x_(X, ii));
% % %
% % %     HesJ =  @(X) HesJ(X) +  Ui*hesarr{1}(x_(X, ii))*(Ui.');
% % % end
% % %

