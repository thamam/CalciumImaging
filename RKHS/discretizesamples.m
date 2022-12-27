function [delta, tkernvec, dsspikesvec, tspikes] = discretizesamples(tvec, spikesvec, deltaTarget)
%DISCRETIZESAMPLES Summary of this function goes here
%   Detailed explanation goes here
% Input - 
%     tvec - uniform time grid on which spikes can be measured
%     spikesvec - indicator vecor same size as tvec if spike occured or not 
%     deltaTarget - this is the discretization resolution we will use for
%     bining the intneisty function estimate (cf. making the Representer
%     theorem hold ...)
% Output - 
%     delta - the actual binning size in case deltaTarget was modified
%     such that delta equals integer number of Ts slots 
%     tkernvec - mx1 vectyor of the time stamps where the kernels will be centered on
%     dsspikesvec - m*1 vector counting the number of spikes in each bin
%     tspikes - the actual time of the spikes

%compute period
Ts = tvec(2) - tvec(1);

%compute downsampling factor
ds = ceil(deltaTarget/Ts);
delta = Ts*ds;

%pad with zeros - no residual
zrpadnum = ds-mod(numel(tvec) , ds);
if zrpadnum>0
    tvec = [tvec(:).', (tvec(end)+Ts*(1:1:zrpadnum))];
    spikesvec = [spikesvec(:).', zeros(1,zrpadnum)];
end

%downsample
tkernvec = delta:delta:tvec(end);
tkernvec=tkernvec(:);
dsspikesvec = sum(reshape(spikesvec,ds,[]),1) ; 
dsspikesvec = dsspikesvec(:) ;

% spikes exact time
tspikes = tvec(logical(spikesvec));

end