function xFit = fit_PoissonGP(t_spikes, x_params, tLims, dt)

% xFit = fit_PoissonGP(t_spikes, x_params, tLims, dt)
%
%
% 2022 - Adam Charles


t_toFit = linspace(min(tLims),max(tLims),ceil((max(tLims)-min(tLims))/dt));

llFunc = @(x) mk_ACJP_ll(x, t_spikes, t_toFit, dt, x_params);

%% Set up and run the Poisson GP fit
options = optimoptions('fminunc','Algorithm','trust-region',...
    'SpecifyObjectiveGradient',true,'HessianFcn','objective');

xFit = fminunc(llFunc, zeros(numel(t_toFit),1), options);

end        
        
function [LL, varargout] = mk_ACJP_ll(x, t_spikes, t_toFit, dt, x_params)

% function LL = mk_ACJP_ll(S, t_sum, n_params, x_params)
%
% Calculate the negative log maginal likelihood value given a particular
% set of GP parameters and PSTH estimate over the summary points.
%
% The cost implemented is
% J(x) = (prod_{t\in t_all} e^{x_t}) * e^{-\Delta \sum_{t\in t_sum} e^{x_t}} * e^{-x^T*k*x/2}
%
% -log(J) = -(1^T*Ker*Krr^{-1}*x_t) + w^Te^{x_t} + x^T*Krr^{-1}*x/2
%
% Grad(-log(J)) = -Krr^{-1}*Ker^T*1 + w.*e^{x_t} + Krr^{-1}*x
%
% Hess(-log(J)) = Delta diag(w.*e^{x_t}) + Krr^{-1}
%
% 2022 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

x    = x(:);               % Ensure that x is a vector
tol  = 1e-2;               % Tolerance for matrix inversion
avWt = 2*ones(size(x));    % Get weights for trapazoidal rule rate integral
avWt([1,end]) = 1;         % Need to set the ends to 1/2
avWt = avWt*dt/2;          % multiply by dt/2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add in PSTH-dependent portions

Kxxb  = mk_GP_mat(t_spikes, t_toFit, x_params);                            % Make the GP cov matrix for the PSTH between the spike times and time points and summary points
Kxbxb = mk_GP_mat(t_toFit, t_toFit, x_params);                             % Make the GP cov matrix for the PSTH between the summary points

Kxxb_sum = sum(Kxxb,1);                                                    % Sum the cross-covariance once (for computational savings)
KIxbxb   = pinv(Kxbxb,tol);                                                % Calculate an inverse once with a given tolerance


LL = 0.5*sum(avWt.*exp(x));                                                % Sum portion of the Poisson term
LL = LL - Kxxb_sum*(KIxbxb*x) + 0.5*x.'*(KIxbxb*x);                        % Add in the contributions depndent on the summary PSTH values

if nargout > 1
    gradLL = -(KIxbxb.')*(Kxxb_sum.') + (KIxbxb*x) +0.5*avWt.*exp(x);  % Gradient has 3 tems
    varargout{1} = gradLL;
end

if nargout > 2
    hessLL = (KIxbxb)  + diag(0.5*avWt.*exp(x));                      % The Hessian only has 2 terms
    varargout{2} = hessLL;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
