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
%% setting simulation and algorithm parameters
% Algorithm parameters
supeps = 1e-3; % support suppers threshold
sig_f = 1; sig_l= 1/2; %kernel initial parameters
k = @(i,j) sig_f*exp(-(i-j).^2./sig_l^2); %kernel handle function
eta = 2; % reg. penalty wieght
gamma = 1e-9; % Tikhonov regularization constant

% Opt. solver params (for Matlab builtin solver)
alg2use = 'quasi-newton';
options = optimoptions('fminunc','Algorithm',alg2use,'SpecifyObjectiveGradient',true);

% options.HessianFcn = 'objective';
% options.CheckGradients = true;
options.FiniteDifferenceType = 'central';
options.FiniteDifferenceStepSize = 1e-10;
options.OptimalityTolerance = 1e-10;
options.MaxIterations = 500;
options.FunctionTolerance = 1e-10;


%% Load data

% change data name into on
% dataname = 'test';
dataname = 'data_082611_cell2_001.mat';
datapath = 'Dataset\CAI_1\GCaMP5k\\processed_data\';

[rawdataout] = loaddata(dataname, datapath) ; %read data for code testing

deltaTarget = 40e-3;

[delta, tdeltavec, dsspikesvec] = discretizesamples(rawdataout.timevec, rawdataout.spikevec, deltaTarget);

%prepare discretized data struct
dataout = rawdataout;
dataout.spikevec = dsspikesvec ;
dataout.timevec = tdeltavec;
dataout.delta = delta;

if 0 %visualize data
    figure(1),clf
    stem(rawdataout.timevec, rawdataout.spikevec, '-*')
    hold all
    stem(tdeltavec,dsspikesvec,'rs')
end

% Compute frame size as 2*tmin, tmin = min_dt s.t. (k(|dt|)<supeps) use
% current sig_l, sig_f, delta check later if need to change frames' size
[tminInd, mintime ] = computeframe(supeps, dataout.delta, k );


% mitime(tminInd) gives the minimum half-time(indices) for local frame size
frmlen = 2*mintime;
m = 2* tminInd;

% cind  - relative index to begining of basisframe  where dataframes
% starts\ends
% m - #basis functions/frame


% Set simulation length as N*frmlen N integer - truncate if needed
Tf = floor(dataout.tf/frmlen); % simulation length

%% Run  streaming solver
[xhat,outstat, XHAT] = runstreamestimate(dataout, m, options, Tf, k, eta, gamma);

saveresults = true;
if saveresults
    save(sprintf('CIResults_%s',nowstamp),'dataout','dataname','datapathg','xhat', 'XHAT', 'options');
end

%% Plot result
% Compute lambdaVec
tplot = dataout.ts:(dataout.delta/100):dataout.tf;
svec = k(tplot(:),glbltvec(:).')*xhat;
lamvec = exp(svec);
%%
figure(1),clf
plot(tplot,lamvec,'--b')
hold all
stem(rawdataout.deltatiks(1:t*m),5*rawdataout.spikevec((1:t*m)))


