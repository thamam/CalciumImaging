function LL = mk_ACJP_ll(t_spikes, x_params, tLims, dt)
% function LL = mk_ACJP_ll(S, t_sum, n_params, x_params)
%
% Calculate the negative log maginal likelihood value given a particular
% set of GP parameters and PSTH estimate over the summary points.
%
% 2022 - Adam Charles


t_toFit = linspace(min(tLims),max(tLims),ceil((max(tLims)-min(tLims))/dt));

llFunc = @(x) mk_ACJP_ll(x, t_spikes, t_toFit, x_params);

%% Set up and run the Poisson GP fit
options = optimoptions('fminunc','Algorithm','trust-region',...
    'SpecifyObjectiveGradient',true,'HessianFcn','objective');

xFit = fminunc(llFunc zeros(numel(t_toFit),1), options);


end        
        
function [LL, varargout] = mk_ACJP_ll(x, t_all, t_sum, x_params)

% function LL = mk_ACJP_ll(S, t_sum, n_params, x_params)
%
% Calculate the negative log maginal likelihood value given a particular
% set of GP parameters and PSTH estimate over the summary points.
%
% 2022 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

x     = x(:);
tol   = 1e-2;
avgWt = 2*ones(size(x));
avgWt([1,end]) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add in PSTH-dependent portions

Kxx   = mk_GP_mat(t_all, t_all, x_params);                                 % Make the GP cov matrix for the PSTH between all the spike time points
Kxxb  = mk_GP_mat(t_all, t_sum, x_params);                                 % Make the GP cov matrix for the PSTH between the spike times and time points and summary points
Kxbxb = mk_GP_mat(t_sum, t_sum, x_params);                                 % Make the GP cov matrix for the PSTH between the summary points

Kxxb_sum = sum(Kxxb,1);                                                    % Sum the cross-covariance once (for computational savings)
KIxbxb   = pinv(Kxbxb,tol);                                                % Calculate an inverse once with a given tolerance

% Cost is
% (sum_{t\in t_all} e^{x_t}) * e^{-\Delta \sum_{t\in t_sum} e^{x_t}} * e^{-x^T*k*x/2}

LL = 0.5*sum(avWt.*exp(x));                                                % Sum portion of the Poisson term
LL = LL - Kxxb_sum*(KIxbxb*x) + 0.5*x.'*(KIxbxb*x);                        % Add in the contributions depndent on the summary PSTH values

if nargout > 1
    gradLL = - (KIxbxb.')*(Kxxb_sum.') + 0.5*0.5*(KIxbxb*x) + 0.5*avWt.*exp(x); % Gradient
    varargout{1} = gradLL;
end

if nargout > 2
    hessLL = 0.5*0.5*(KIxbxb)  + diag(0.5*avWt.*exp(x));       % Hessian
    varargout{2} = gradLL;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
