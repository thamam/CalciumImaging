%% Streaming spikes regression - main frame
%mainframe.m
% written by Tomer H.Hamam tomerhamam@gmail.com 09.16.2021
% streaming estimation of intensity function from spikes measurements
% User options:
% Select data to process by setting the following parameters:
% dataname = 'test' / 'user' /
% datapath
%
%%%%%
close all;
clear;
clc;
addpath(genpath('DataGen/'));
addpath(genpath('RKHS'));
addpath(genpath('Types'));

saveresults = true;
makevideo = true;
%% setting simulation and algorithm parameters
%% BASIC USER INTERFACE

datatype  = DatasetsType.Sim;    

% Select data: % change data name into dataname = 'test' , or name of of
% the CI datasets files, and then set the path to dat in 'datapath'
if isequal(datatype, DatasetsType.GCaMP5k)
    data_options.dataname = 'data_080511_cell7_002.mat';
    data_options.datapath = 'Datasets\GCaMP5k\\processed_data\';
elseif isequal(datatype, DatasetsType.Sim)
    data_options = createProblemStruct();    
    data_options.tmax = 12;
end
 
BUFFERLENGTH = 2;

% prepare animated video of results
MAKEVIDEO = false;

% Start interactive plot of rsults
SMARTPLOT = true;

% Option to limit length of simulation
maxsimleng = 20; %inf; % in #of frames

% frame size scale factor (default m = 2*tmin |tmin = min_delt
% k(|t-delt|)<supeps
fmag = 1;
%% ADVANCED USER INTERFACE


% kernel parameters
supeps = 1e-3; % support suppers threshold

sig_f = 1; sig_l= 1/2; %kernel initial parameters

eta = 2; % penalty wieght
gamma = 1e-9; % Tikhonov regularization constant

% kernels spacings (and accuracy of computing integral in ML )
% maybe rounded to align with time resolution of avaiable data
deltaTarget = 40e-3;

% Choose which algorithm to use: 'trust-region' , or  'quasi-newton'
alg2use = 'trust-region' ;
% alg2use = 'quasi-newton';

% Plotting resolution
TsPlotfactor = 10; % TsPlot  = delta/TsPlotfactor;

% Can set to false and provide mat file that has 'xhat_batch' already saved
recomputebatch = true; 
% xbatchfilename  =[];



%%
% Opt. solver params (for Matlab builtin solver)
% Algorithm parameters
options = optimoptions('fminunc','Algorithm',alg2use,'SpecifyObjectiveGradient',true);
options.HessianFcn = 'objective';
% options.CheckGradients = true;
% options.FiniteDifferenceType = 'central';
% options.FiniteDifferenceStepSize = 1e-10;
% options.OptimalityTolerance = 1e-10;
% options.MaxIterations = 500;
% options.FunctionTolerance = 1e-10;

%% Load data
[rawdataout] = loaddata(datatype, data_options) ; %read data for code testing
[delta, tkernvec, dsspikesvec, tspikes] =...
    discretizesamples(rawdataout.timevec, rawdataout.spikevec, deltaTarget);

%prepare discretized data struct
dataout = rawdataout;
dataout.binsspikes = dsspikesvec ; %spikes counted per dicrete bin
dataout.timevec = tkernvec; % discretization time marks
dataout.delta = delta;
dataout.tspikes = tspikes;

if 0 %visualize data
    figure(1),clf
    stem(rawdataout.timevec, rawdataout.spikevec, '-*')
    hold all
    stem(tkernvec,dsspikesvec(1:end-1),'rs')
end

%Define kernel handle function - square exponential 
k = @(i,j) sig_f*exp(-(i-j).^2./sig_l^2); 

% Compute frame size as 2*tmin, tmin = min_dt s.t. (k(|dt|)<supeps) use
% current sig_l, sig_f, delta check later if need to change frames' size
[tminInd, mintime ] = computeframe(supeps, dataout.delta, k );


% Determine frame length : mintime(tminInd) gives the minimum
% half-time(indices) for local frame size
frmlen = fmag*2*mintime;
m = fmag*2*tminInd;

% cind  - relative index to begining of basisframe  where dataframes
% starts\ends
% m - #basis functions/frame


% Set simulation length as N*frmlen N integer - truncate if needed
MF = floor(dataout.tf/frmlen); % simulation length in frames

if MF > maxsimleng
    MF = floor(maxsimleng);
end

%% Run  streaming solver
[xhat,outstat,ktvec, dtvec, XHAT, tauArray]= ...
    runstreamestimate(dataout, m, options, MF, k, eta, gamma, BUFFERLENGTH);

nowstamp = datetime('now','TimeZone','local','Format','d-MMM-y_HH_mm');
TsPlot  = delta/TsPlotfactor;

%% bACK Populate XHAt solutions 
if BUFFERLENGTH>MF
    XHAT_full = XHAT;
else
    XHAT_full=XHAT;
    for ii=BUFFERLENGTH+1:MF
        missingInd = 1:(ii-BUFFERLENGTH)*m;
        XHAT_full(missingInd,ii) = XHAT_full(missingInd,ii-1);
    end
end


%% Manual plot  - for debug
plotxhat_tf=1;
if plotxhat_tf
    xhat = XHAT_full(:,end);
    %% Compute lambdaVec
    tplot = dataout.ts:(dataout.delta/100):ktvec(end);
    lamvec = reconstlambda(k, xhat,ktvec,tplot);
    figure(1),clf
    plot(tplot,lamvec,'b','LineWidth',1.5)
    hold all
    lambatspikes = reconstlambda(k, xhat,ktvec,dataout.tspikes)
    stem(dataout.tspikes,lambatspikes,'ks','MarkerSize',5,'LineWidth',0.5,'LineStyle',':')
end

%% Compute batch(off-line) solve
if recomputebatch
    X0 = XHAT_full(:,end);
    [xstar,fval, exitflag] = globalsolver(X0, dataout, ktvec, k,eta);
else
    xstar = load(xbatchfilename,'xhat_batch');
end

if saveresults
    save(sprintf('CIResults_%s',nowstamp),'dataout','datatype',...
        'data_options','xhat', 'XHAT', 'options','k','MF','ktvec', 'TsPlot','xstar');
end

%% Prepare results video animation file
if (MAKEVIDEO)
    vidtit = sprintf('Streamspikes_%s_%s',dataname, nowstamp);
    makespikesvideo(vidtit, XHAT_full, k, m,  MF, ktvec, TsPlot, tauArray, xstar )    ;
end
%%
%% Smart pllot
% load('CIResults_20-Sep-2021_16_33.mat');

%% Load missing variables for playing saved results
if SMARTPLOT
    smartplot(XHAT_full, k, m,  MF, ktvec, TsPlot, tauArray, BUFFERLENGTH, xstar);
end


