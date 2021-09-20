function [dataout] = loaddata(dataname, datapath)
%loaddata Load spikes data
%   Detailed explanation goes here


switch dataname
    case 'test'
        Ts = 40e-3;
        [ts, tf] = deal(0, 15);
        timevec = ts+Ts:Ts:tf;
        spikevec = randi([0,1],numel(timevec),1);
        dataout.spikevec = spikevec(:);
        dataout.Ts=Ts;
        dataout.ts=ts;
        dataout.tf=tf;
        dataout.timevec = timevec(:);        
    otherwise
        obj_ = load(strcat(datapath,dataname),'obj');
        obj=obj_.obj;
        % Fluorocensce data - ignore for now
        %fmean_roi=obj.timeSeriesArrayHash.value{1}.valueMatrix;
        %fmean_neuropil=obj.timeSeriesArrayHash.value{2}.valueMatrix;
        %fmean_comp=fmean_roi-0.7*fmean_neuropil;        
%         t_frame=obj.timeSeriesArrayHash.value{1}.time;        
%         filt=obj.timeSeriesArrayHash.value{4}.valueMatrix;
        t_ephys=obj.timeSeriesArrayHash.value{4}.time;       
        
        spikevec=obj.timeSeriesArrayHash.value{5}.valueMatrix;                      
        dataout.spikevec = spikevec(:);
        dataout.Ts = t_ephys(2)- t_ephys(1);
        dataout.ts = t_ephys(1) ;
        dataout.tf = t_ephys(end);
        dataout.timevec = t_ephys(:);    
end


end
