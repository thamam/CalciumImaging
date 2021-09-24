function []=smartplot(XHAT, k, m,  MF, ktvec, TsPlot, tauArray, BUFFERLENGTH)
%SMARTPLOT Summary of this function goes here
%   Detailed explanation goes here
% tspikes = dataout.tspikes;
resnum = size(XHAT,2);
xstar = XHAT(1:(MF*m),end);
tvecplotstar = (ktvec(1):TsPlot:ktvec(end)).';
lambdastar = reconstlambda(k, xstar, ktvec, tvecplotstar);
% Generate a set of frames, get the frame from the figure, and then write each frame to the file.
t=1;
dirkey = 'u';
    h1 = figure(1);
    clf(h1);
    hup = subplot(2,1,1);
    hold all

    hdown = subplot(2,1,2);
    hold all
    
while(~strcmp(dirkey,'e'))   
    figure(h1)
    %Plotting solution up to time t in resnum
    xhat_t = XHAT(1:t*m,t);
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
        %Print global scale
        plot(hup, tvecplotstar, lambdastar,'LineStyle',':','Color',[0.1,0.1,0.1],'LineWidth',1);
        stem(hup, tau_t, lamfunc(tau_t.'),'ks')
        plot(hup, tvecplot, lamhat,'LineWidth',2);
        
        %Print local scale
%         locind = 1:m;
%         loctvec = ktvec_t(locind);
%         localTvec = ktvec_t(1):TsPlot:ktvec_t(end);
%         loclamhat = reconstlambda(k, xhat_t, ktvec_t, tvecplot)
        
        stem(hdown, tau_t, lamfunc(tau_t.'),'ks')
        plot(hdown, tvecplot, lamhat,'LineWidth',2);
        
%         stem(hdown, tau_t, lamfunc(tau_t.'),'ks')
%         plot(hdown, loctvec, loclamhat,'LineWidth',2);        
        hdown.YLim = [0 max(lamhat)*1.05];
        hdown.XLim= [tvecplot(1) tvecplot(end)];
    else
        %Print global scale
        hold(hup,'off')
        plot(hup, tvecplotstar, lambdastar,'LineStyle',':','Color',[0.1,0.1,0.1],'LineWidth',1);        
        hold(hup,'on')
        %plotting previous estimate
        plot(hup, tvecplot_c,lam_c,'g', 'linewidth',1,'LineStyle','--')        
        %plotting new estiamte
        stem(hup, tau_t, lamfunc(tau_t.'),'ks','filled','MarkerFaceColor','y');
        plot(hup, tvecplot, lamhat,'b','LineWidth',1.5);
        
        %Print local scale
        hold(hdown,'off')

         %plotting previous estimate
        plot(hdown, tvecplot_c,lam_c,'g', 'linewidth',1,'LineStyle','--')        
        %plotting new estiamte
        stem(hdown, tau_t, lamfunc(tau_t.'),'ks','filled','MarkerFaceColor','y');
        plot(hdown, tvecplot, lamhat,'b','LineWidth',1.5);
        
        roiInd = [ max(0,tvecplot(end)-(BUFFERLENGTH)*(ktvec(2)-ktvec(1))*m) tvecplot(end)];        
        lamroiind = (tvecplot>roiInd(1)&tvecplot<roiInd(2));
        hdown.YLim = [0 max(lamhat(lamroiind))*1.15];
        hdown.XLim= roiInd;      
    end
    %save current estimate for future plotting
    lam_c = lamhat;
    tvecplot_c = tvecplot;
    
    %% figure anotation
    title(hup,sprintf('Streaming estimate: $$ {\\lambda}(t), t_f = %.2f~[s]$$, \\quad  T(frame)= %i',...
        ktvec_t(end-m/2), t),'FontSize',18,'FontName','Sitka Small','Interpreter','latex');        
    xlabel('time[t]','FontSize',20)
    ylabel('$$\hat{\lambda} (t) $$','Interpreter','latex','FontSize',20)
    
    if t==MF
        title(hdown,sprintf('T(frame \\#) = %i  EOF', t),'FontSize',18,...
            'FontName','Sitka Small','Interpreter','latex'); 
    else
        title(hdown,sprintf('T(frame \\#) = %i', t),'FontSize',18,...
            'FontName','Sitka Small','Interpreter','latex'); 
    end
    %% manage user keyboard interface
    dirkey=input('Press f (Fwd) / b (Bwd) / e (Esc)\n','s');
    if (strcmp(dirkey, 'f'))
        t = min(t+1,MF);
    elseif(strcmp(dirkey, 'b'))
        t = max(t-1,1);
    end
    
   
       
    
    
end


end

