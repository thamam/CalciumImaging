function K = mk_GP_mat(T1, T2, params)

% K = mk_GP_mat(T1, T2, params)
% 
% Function to make a GP covariance matrix. 
% 
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if numel(params) == 1
    params = [1, params, 2];
elseif numel(params) == 2
    params = [params(1), params(2), 2];
elseif numel(params) == 3
    params = params(:).';
else
    warning('params has too many elements. Using the first 3 only!')
    params = [params(1), params(2), params(3)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate matrix

K = abs(bsxfun(@plus, T1(:), -T2(:).'));                                   % Get differences b/w each pait of time points
K = params(1)*exp(-(K.^(params(3)))./(params(2).^params(3)));              % calculate rho*e^{-|dt|^p/l^p}

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%