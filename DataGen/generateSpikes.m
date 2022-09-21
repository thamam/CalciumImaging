function datOut = generateSpikes(varargin)

% datOut = generateSpikes(varargin)
% 
% Function to generate spike events from a latent GP rate. 
%
% Example:
%
% dataOut = generateSpikes('tmax',6,'N_trial',5,'rateOffset',1)
%
% 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up simulation parameters

p = inputParser; % create input scheme
addParameter(p, 'tmax'      ,     4             );                         % Time window length (in seconds) to simulate spikes over
addParameter(p, 'x_params'  ,     [8, 1, 2]     );                         % Parameters of the latent rage Gaussian Process
addParameter(p, 'N_samp'    ,     50            );                         % Samples are in per-second units
addParameter(p, 'N_trial'   ,     10            );                         % Number of trials to simulate (resampling of spikes given the same rate)
addParameter(p, 'rateOffset',     0             );                         % Possible offset for the rage (if a nonzero mean is required)
% Parameters below are not recommended to change
addParameter(p, 'cosntx'    ,     false         );
addParameter(p, 'overDisp'  ,     false         ); % normalize option
addParameter(p, 'threshVal' ,     0             );
addParameter(p, 'min_dt'    ,     2e-3          );
addParameter(p, 'pmax'      ,     Inf           );
addParameter(p, 'plotOpt'   ,     false         );
addParameter(p, 'n_params'  ,     [0.5, 0.5, 1] );
addParameter(p, 'sampOpt'   ,     'circest'     );

parse(p,varargin{:});                                                      % Parse and validate input arguments(contained in cell array varargin)
p = p.Results;                                                             % Contains the validated values of the inputs.

p.N_samp = p.tmax*p.N_samp;                                                % Samples are in per-second units
tt_samp  = linspace(0,p.tmax,p.N_samp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a set of trials using a single PSTH

if p.cosntx
    x_true = sqrt(x_params(1))*randn(1);
    x_proj = x_true/x_params(1);
else
    Kxx_full = mk_GP_mat(tt_samp, tt_samp, p.x_params);                    % Make GP matrices for the data
    x_true   = Inf;
    while (max(exp(x_true)) > 5)||(max(exp(x_true)) < exp(1))
        switch lower(p.sampOpt)
        case 'circest'
            x_true = sampRandMVNGP(zeros(numel(tt_samp),1),Kxx_full(1,:).');
        otherwise
            x_true = mvnrnd(zeros(numel(tt_samp),1),Kxx_full);             % Resample until the max values are not insane
        end
    end 
    x_true = x_true + p.rateOffset;
    x_proj = pinv(Kxx_full,1e-4)*x_true(:);
end

%%

n_true = cell(p.N_trial,1);
n_proj = cell(p.N_trial,1);
p_rate = cell(p.N_trial,1);
evt    = cell(p.N_trial,1);

if p.overDisp
    Knn_full = mk_GP_mat(tt_samp, tt_samp, p.n_params);                    % Make GP matrices for the noise
    for kk = 1:p.N_trial
        n_true{kk} = mvnrnd(-0.5*p.n_params(1)*ones(numel(tt_samp),1),...
                                                                 Knn_full);
        if p.cosntx
            rmax = 1*max(exp(x_true + n_true{kk}(:)));
            n_proj{kk} = pinv(Knn_full,1e-4)*(n_true{kk}(:) + ... 
                              0.5*p.n_params(1)*ones(size(n_true{kk}(:))));
            p_rate{kk} = @(t) min(exp(-0.5*p.n_params(1) ... 
                         + mk_GP_mat(t, tt_samp, p.n_params)*n_proj{kk} ... 
                                                          + x_true), rmax);
        else
            rmax = 1*max(exp(x_true(:) + n_true{kk}(:)));
            n_proj{kk} = pinv(Knn_full,1e-4)*(n_true{kk}(:) + ...
                              0.5*p.n_params(1)*ones(size(n_true{kk}(:))));
            p_rate{kk} = @(t) min(exp(-0.5*p.n_params(1) ...
                        + mk_GP_mat(t, tt_samp, p.n_params)*n_proj{kk} ...
                        + mk_GP_mat(t, tt_samp, p.x_params)*x_proj), rmax);
        end
    end
else
    for kk = 1:p.N_trial
        n_true{kk} = zeros(size(x_true));
        if p.cosntx
            rmax = 1*max(exp(x_true));
            p_rate{kk} = @(t) min(exp(x_true), rmax);
        else
            rmax = 1*max(exp(x_true(:)));
            p_rate{kk} = @(t) min(exp(...
                          mk_GP_mat(t, tt_samp, p.x_params)*x_proj), rmax);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample continuous time spikes using Ogata

for kk = 1:p.N_trial
    min_dt     = 0;
    while min_dt < p.min_dt
        evt{kk} = ogatapoisson(p_rate{kk},rmax,p.tmax,p.pmax);
        if numel(evt{kk}) > 1
            min_dt = min(diff(evt{kk}));
        else
            min_dt = Inf;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up output

datOut.x_true = x_true;
datOut.n_true = n_true;
datOut.evt    = evt;
datOut.p_rate = p_rate;
datOut.rmax   = rmax;
datOut.min_dt = min_dt;
datOut.mean_proc = mean(exp(bsxfun(@plus, x_true(:), cell2mat(n_true).')),2);
datOut.std_proc  = std(exp(bsxfun(@plus, x_true(:), cell2mat(n_true).')),[],2);
datOut.tt_samp = tt_samp; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optional plotting

if p.plotOpt
    trial_sel = randsample(p.N_trial,min(p.N_trial,5));
    subplot(1,2,1), cla
    subplot(1,2,1), hold on
    plot(tt_samp, exp(x_true(:).*ones(numel(tt_samp),1)), 'LineWidth', 3)
    for kk = 1:numel(trial_sel)
        plot(tt_samp, exp(n_true{trial_sel(kk)}(:)), 'k', 'LineWidth', 3)
        plot(tt_samp, exp(x_true(:) + n_true{trial_sel(kk)}(:)), 'r', 'LineWidth', 3)
    end
    ylabel('Poisson rate', 'FontSize', 20)
    box off
    legend('PSTH', 'Dispersion noise', 'Process rate')
    set(gca, 'FontSize', 18, 'TickDir', 'out', 'XLim', [0,p.tmax])
    subplot(1,2,1), hold off


    mean_proc = mean(exp(bsxfun(@plus, x_true(:), cell2mat(n_true).')),2);
    std_proc  = std(exp(bsxfun(@plus, x_true(:), cell2mat(n_true).')),[],2);
    xlabel('Time (s)', 'FontSize', 20)
    subplot(1,2,2), cla
    subplot(1,2,2), hold on
    plot(tt_samp, [exp(x_true(:).*ones(numel(tt_samp),1)), mean_proc], 'LineWidth', 3)
    plot(tt_samp, [mean_proc-std_proc, mean_proc+std_proc], 'r', 'LineWidth', 3)
    for kk = 1:numel(trial_sel)
    %     plot(tt_samp, [exp(n_true{trial_sel(kk)}(:)), exp(x_true(:) + n_true{trial_sel(kk)}(:))], 'r', 'LineWidth', 2)
        stem(evt{trial_sel(kk)}, mean(mean_proc)*ones(size(evt{trial_sel(kk)})), 'LineWidth', 3)
    end
    % stem(evt{trial_sel}, mean(mean_proc)*ones(size(evt{trial_sel})), '.-k', 'LineWidth', 3)
    subplot(1,2,2), hold off
    xlabel('Time (s)', 'FontSize', 20)
    box off
    set(gca, 'FontSize', 18, 'TickDir', 'out', 'XLim', [0,p.tmax], 'YTick', [])

    set(gcf, 'color', [1,1,1])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x_true = sampRandMVNGP(u,sigVec)

N = numel(sigVec);

sigVec = sigVec(:)';
sigVec = [sigVec,zeros(1,N),fliplr(sigVec(2:end))]./numel(sigVec);


x_true = randn(size(u)); % Initial simulation of iid data
x_true = sqrt(fft(sigVec(:))).*[x_true(:);zeros(N,1);zeros(N-1,1)];
x_true = real(fft((x_true)));
x_true = x_true(1:N);
x_true = x_true + u;     % Offset by the mean

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
