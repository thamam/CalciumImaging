function x_out = GPbaseline(tSpikes, tmax, nPts, sig, varargin)


% Generate a Gaussian
tt = linspace(0,tmax,nPts);

GaussFun = @(x) exp( (tt-x).^2/(sig^2) )

%% Convolve with the spike-trains

x_out = 0;
for ii = 1:numel(tSpikes)
    x_out = x_out + GaussFun(t_Spikes);
end

end
