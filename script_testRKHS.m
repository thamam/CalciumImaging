%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function to test the RKHS solver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up simulation data
% 'tmax'     , 4             
% 'cosntx'   , false         
% 'overDisp' , false         
% 'threshVal', 0             
% 'N_samp'   , 50            
% 'N_trial'  , 10            
% 'pmax'     , Inf           
% 'plotOpt'  , false         
% 'x_params' , [8, 1, 2]     
% 'min_dt'   , 2e-3 
pars.tmax     = 6;
pars.x_params = [8, 1, 2]; % [p1,p2,p3] GP Kernel variables K = p1 *exp(|| dt||^p3/(p2))

datOut = generateSpikes('N_trial', 10, 'tmax', pars.tmax, 'x_params', pars.x_params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Poisson-GP (offline)

xFit = fit_PoissonGP(cell2mat(datOut.evt), pars.x_params , [0,pars.tmax], 1e-2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Online inference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot results

figure(1), cla
figure(1), hold on
% plot(linspace(0,pars.tmax,numel(datOut.x_true)), datOut.x_true, '-b',...
%     linspace(0,pars.tmax,numel(xFit)), xFit, '-r')
plot(linspace(0,pars.tmax,numel(datOut.x_true)), exp(datOut.x_true), '--b',...
    linspace(0,pars.tmax,numel(xFit)), exp(xFit)/20, '--r')
stem(cell2mat(datOut.evt), ones(size(cell2mat(datOut.evt))))
figure(1), hold off