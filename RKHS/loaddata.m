function [dataout] = loaddata(datatype, options)
%loaddata Load spikes data
%   Detailed explanation goes here

switch datatype
    case DatasetsType.Synt
        Ts = 40e-3;
        [ts, tf] = deal(0, 15);
        timevec = ts+Ts:Ts:tf;
        spikevec = randi([0,1],numel(timevec),1);
        dataout.spikevec = spikevec(:);
        dataout.Ts=Ts;
        dataout.ts=dataout.spikevec(1);
        dataout.tf=dataout.spikevec(end);
        dataout.timevec = timevec(:);
        dataout.src = 'sim';

    case DatasetsType.Sim
        t_final = options.tmax;
        t_start = 0;
        dataOut = generateSpikes('tmax',options.tmax,'N_trial', options.N_trial,'rateOffset', options.rateOffset);
        dataout.ts = t_start;
        dataout.tf = t_final;
        dataout.tspikes = sort(vertcat(dataOut.evt{:}),'ascend');                       
        dt = (dataOut.tt_samp(2)- dataOut.tt_samp(1))/10;
        tt_vec = linspace(t_start,t_final,ceil((t_final-t_start)/dt)).';
        [~, index]  = min(abs(tt_vec-dataout.tspikes.'));
        dataout.spikevec = zeros(size(tt_vec));
        dataout.spikevec(index) = 1;
        dataout.timevec = tt_vec;
        dataout.Ts = tt_vec(2)- tt_vec(1);  % TODO (Tomer) do we need it?
        dataout.src = 'sim';

    case DatasetsType.GCaMP5k
        obj_ = load(strcat(options.datapath ,options.dataname ),'obj');
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
        dataout.src = 'GC';

end


end
