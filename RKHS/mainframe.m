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

saveresults = true;
makevideo = true;
%% setting simulation and algorithm parameters
% Select data
% change data name into on

% dataname = 'test';
% datapath=[];

% addpath(genpath('C:\Users\tomer\MATLAB\Projects\CI\Datasets'));

dataname = 'data_080511_cell7_002.mat';
datapath = '..\Datasets\GCaMP5k\\processed_data\';

%%
% Algorithm parameters
supeps = 1e-3; % support suppers threshold
sig_f = 1; sig_l= 1/2; %kernel initial parameters
k = @(i,j) sig_f*exp(-(i-j).^2./sig_l^2); %kernel handle function
eta = 2; % reg. penalty wieght
gamma = 1e-6; % Tikhonov regularization constant

% Opt. solver params (for Matlab builtin solver)
% alg2use = 'quasi-newton';
alg2use = 'trust-region' ;
options = optimoptions('fminunc','Algorithm',alg2use,'SpecifyObjectiveGradient',true);

options.HessianFcn = 'objective';
% options.CheckGradients = true;
% options.FiniteDifferenceType = 'central';
% options.FiniteDifferenceStepSize = 1e-10;
% options.OptimalityTolerance = 1e-10;
% options.MaxIterations = 500;
% options.FunctionTolerance = 1e-10;

%% Load data
[rawdataout] = loaddata(dataname, datapath) ; %read data for code testing
deltaTarget = 40e-3;
[delta, tkernvec, dsspikesvec, tspikes] = discretizesamples(rawdataout.timevec, rawdataout.spikevec, deltaTarget);

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
    stem(tkernvec,dsspikesvec,'rs')
end

% Compute frame size as 2*tmin, tmin = min_dt s.t. (k(|dt|)<supeps) use
% current sig_l, sig_f, delta check later if need to change frames' size
[tminInd, mintime ] = computeframe(supeps, dataout.delta, k );


% mitime(tminInd) gives the minimum half-time(indices) for local frame size
fmag = 1;
frmlen = fmag*2*mintime;
m = fmag*2*tminInd;

% cind  - relative index to begining of basisframe  where dataframes
% starts\ends
% m - #basis functions/frame


% Set simulation length as N*frmlen N integer - truncate if needed
MF = floor(dataout.tf/frmlen); % simulation length in frames

%% Run  streaming solver
 [xhat,outstat,ktvec, dtvec, XHAT, tauArray]= runstreamestimate(dataout, m, options, MF, k, eta, gamma);

nowstamp = datetime('now','TimeZone','local','Format','d-MMM-y_HH_mm');
if saveresults
    save(sprintf('CIResults_%s',nowstamp),'dataout','dataname','datapath','xhat', 'XHAT', 'options');
end



%% Plot result
plotxhat_tf=0;
if plotxhat_tf
    %% Compute lambdaVec
    tplot = dataout.ts:(dataout.delta/100):ktvec(end);
    lamvec = reconstlambda(k, xhat,ktvec,tplot);    
    figure(1),clf
    plot(tplot,lamvec,'--b')
    hold all
    stem(dataout.timevec(1:MF*m),5*dataout.spikevec((1:MF*m)))
end
if (makevideo == true)
    TsPlot = delta/10;
    vidtit = sprintf('Streamspikes_%s_%s',dataname, nowstamp);
    makespikesvideo(vidtit, XHAT, k, m,  MF, ktvec, TsPlot, tauArray)    ;
end
%%
%% Smart pllot
load('CIResults_20-Sep-2021_16_33.mat');
smartplot(XHAT, k, m,  MF, ktvec, TsPlot, tauArray);
