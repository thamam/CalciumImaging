function []=smartplot(XHAT, k, m,  MF, ktvec, TsPlot, tauArray)
%SMARTPLOT Summary of this function goes here
%   Detailed explanation goes here
% tspikes = dataout.tspikes;


resnum = size(XHAT,2);
xstar = XHAT{resnum};
tvecplotstar = (ktvec(1):TsPlot:ktvec(end)).';
lambdastar = reconstlambda(k, xstar, ktvec, tvecplotstar);
% Generate a set of frames, get the frame from the figure, and then write each frame to the file.

t=1;
dirkey = 'u';
while(~strcmp(dirkey,'e'))
    h1 = figure(1);
    clf(h1);
    hold all
    
    
    %Plotting solution up to time t in resnum
    xhat_t = XHAT{t};
    Nt = numel(xhat_t);
    ktvec_t = ktvec(1:Nt);
    tvecplot = (ktvec_t(1):TsPlot:ktvec_t(end));
    tvecplot = tvecplot(:);
    tau_t = tauArray{t};
    tau_t = tau_t(:);
    
    lamfunc = @(t) reconstlambda(k, xhat_t, ktvec_t, t);
    lamhat = reconstlambda(k, xhat_t, ktvec_t, tvecplot);
    % Plotting solution estiamtes
    
    if t==1
        plot(tvecplotstar, lambdastar,'LineStyle',':','Color',[0.1,0.1,0.1],'LineWidth',1);
        stem(tau_t, lamfunc(tau_t.'),'ks')
        plot(tvecplot, lamhat,'LineWidth',2);
    else
        clf(h1);
        plot(tvecplotstar, lambdastar,'LineStyle',':','Color',[0.1,0.1,0.1],'LineWidth',1);
        
        hold all;
        %plotting previous estimate
        plot(tvecplot_c,lam_c,'k', 'linewidth',1,'LineStyle','--')
        
        %plotting new estiamte
        stem(tau_t, lamfunc(tau_t.'),'ks','filled','MarkerFaceColor','y');
        plot(tvecplot, lamhat,'LineWidth',2);
    end
    
    %save current estimate
    lam_c = lamhat;
    tvecplot_c = tvecplot;
    
    %     axis square
    
    %     h1.Children.XLim = [-12,12];
    %     h1.Children.YLim = [-12,12];
    %%
    title(sprintf('Rate function streaming reconstruction: $$ t_f = %.2f ~ [s]$$',...
        ktvec_t(end-m/2)),'FontSize',18,'FontName','Sitka Small','Interpreter','latex');
    
    
    %%
    xlabel('time[t]','FontSize',20)
    ylabel('$$\hat{\lambda} (t) $$','Interpreter','latex','FontSize',20)
    
    
    dirkey=input('Press f (Fwd) / b (Bwd) / e (Esc)\n','s');
    if (strcmp(dirkey, 'f'))
        t = min(t+1,MF);
    elseif(strcmp(dirkey, 'b'))
        t = min(t-1,1);
    end
end


end

