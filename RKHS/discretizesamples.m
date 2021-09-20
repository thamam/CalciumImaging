function [delta, tkernvec, dsspikesvec, tspikes] = discretizesamples(tvec, spikesvec, deltaTarget)
%DISCRETIZESAMPLES Summary of this function goes here
%   Detailed explanation goes here


%compute period
Ts = tvec(2) - tvec(1);

%compute downsampling factor
ds = ceil(deltaTarget/Ts);
delta = Ts*ds;

%pad with zeros - no residual
zrpadnum = mod(numel(tvec) , ds);
if zrpadnum>0
    tvec = [tvec(:).', zeros(1,zrpadnum)];
    spikesvec = [zrpadnum(:).', zeros(1,zrpadnum)];
end

%downsample
tkernvec = delta:delta:tvec(end);
tkernvec=tkernvec(:);
dsspikesvec = sum(reshape(spikesvec,ds,[]),1) ; 
dsspikesvec = dsspikesvec(:) ;

% spikes exact time
tspikes = tvec(logical(spikesvec));

end