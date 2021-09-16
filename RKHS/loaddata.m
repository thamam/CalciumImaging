function [dataout] = loaddata(dataname, datapath)
%loaddata Load spikes data 
%   Detailed explanation goes here


switch dataname
    case 'test'
        delta = 40e-3;
        [ts, tf] = deal(0, 15);
        deltatiks = ts+delta:delta:tf;
        spikevec = randi([0,1],numel(deltatiks),1);
        dataout.spikevec = spikevec;
        dataout.delta=delta;
        dataout.ts=ts;
        dataout.tf=tf;
        dataout.deltatiks = deltatiks;
    case 'user'
        
        
end


end




