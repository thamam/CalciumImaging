%% Example code to read cell-attached recording file format
%
%  this .m file contains example code to read .mat file containing
%  GCaMP5,6s,6f in vivo imaging/ephys data published in 
%  (Chen et. al. 2013 Nature; Akerboom, Chen 2012 J. Neurosci)
%
% Ephys data were recorded at 10 KHz
% Imaging data were recorded at 60 Hz % ????time indices shows 50 Hz
%
% Each .mat data file contains a variable named 'obj'
% key recording traces and time base can be accessed by the following:
%
% traces=obj.timeSeriesArrayHash.value{id}.valueMatrix
% time=  obj.timeSeriesArrayHash.value{id}.time
%
% id=
% 
% 1: fmean_roi
% 2: fmean_neuropil
% 3: raw_ephys
% 4: filtered_ephys
% 5: detected_spikes
%
%
%  Tsai-Wen Chen, 2015/01/27

%%
% file=uigetfile('data*.mat');
% load(file);
% change to read the recording of choice
% load('data_071411_cell1_005.mat');
% load('data.mat');
% load('data_080311_cell2_001.mat'); % clean but very small number of spikes
% load('data_080511_cell7_002.mat');
% load('data_080511_cell12_002.mat');
% load('data_082611_cell1_002.mat'); % clear diagram - second choice
load('data_082611_cell2_001.mat');% clear diagram - first choice
% load('data_090111_cell1002.mat'); % clean with some outliers or falsealarms
% load('data_090711_cell4003.mat');
% load('data.mat');


%---------------------------------------

fmean_roi=obj.timeSeriesArrayHash.value{1}.valueMatrix;
fmean_neuropil=obj.timeSeriesArrayHash.value{2}.valueMatrix;
fmean_comp=fmean_roi-0.7*fmean_neuropil;
t_frame=obj.timeSeriesArrayHash.value{1}.time;

filt=obj.timeSeriesArrayHash.value{4}.valueMatrix;
t_ephys=obj.timeSeriesArrayHash.value{4}.time;

detected_spikes=obj.timeSeriesArrayHash.value{5}.valueMatrix;
spike_time=t_ephys(detected_spikes);

%%
figure;
h1=subplot(2,1,1);
plot(t_frame,fmean_comp)

h2=subplot(2,1,2);
plot(t_ephys,filt,'r')
hold on;
plot(spike_time,filt(detected_spikes),'.k');
linkaxes([h1,h2],'x');
title('Data visualization')


