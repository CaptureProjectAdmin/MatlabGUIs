% ITPCSimTest_run
% Simulation: theta (4-8 Hz) sine bursts with phase jitter, ITPC + power
% computed and plotted with the same method/layout as MultTransITPCEyeGUI /
% RWAnalysis2.plotMultTransITPCEye.
%
% Does not modify any existing codebase files. Requires Wavelet Toolbox
% (morseSpecGram / cwtfilterbank) and Statistics Toolbox if CorrectionType='fdr'.
%
% Usage:
%   >> ITPCSimTest_run

clear; close all; clc;

%% Parameters (edit freely)
cfg = struct();
cfg.Fs = 250;                    % Hz (Neuropace-like)
cfg.DurationSec = 6;             % total window
cfg.OnsetSec = 4;                % burst appears ~second 3
cfg.NTrials = 50;                % number of sine-wave trials
cfg.FreqHz = 6;                  % carrier (or [4 8] for random theta)
cfg.PhaseJitterRad = pi;       % +/- 60 deg phase jitter across trials
cfg.BurstHalfWidthSec = 0.75;    % Hann half-width around onset
cfg.Amplitude = 30;              % signal amplitude
cfg.NoiseStd = 10;               % noise std
cfg.Seed = 1;

cfg.NPerm = 200;                 % raise toward 500-1000 for publication-like nulls
cfg.PVal = 0.05;
cfg.CorrectionType = 'pixel';    % 'pixel' or 'fdr'
cfg.CLim = [-10 10];

%% Generate jittered theta trials
[sig, t, meta] = ITPCSimTest_generateSignals( ...
    'Fs',cfg.Fs, ...
    'DurationSec',cfg.DurationSec, ...
    'OnsetSec',cfg.OnsetSec, ...
    'NTrials',cfg.NTrials, ...
    'FreqHz',cfg.FreqHz, ...
    'PhaseJitterRad',cfg.PhaseJitterRad, ...
    'BurstHalfWidthSec',cfg.BurstHalfWidthSec, ...
    'Amplitude',cfg.Amplitude, ...
    'NoiseStd',cfg.NoiseStd, ...
    'Seed',cfg.Seed);

fprintf('Generated %d trials | mean f=%.2f Hz | phase jitter +/-%.0f deg\n', ...
    meta.NTrials, mean(meta.FreqHz), rad2deg(meta.PhaseJitterRad));

%% Compute ITPC + power (same pipeline as RWAnalysis2)
S = ITPCSimTest_compute(sig, cfg.Fs, ...
    't', t, ...
    'NPerm', cfg.NPerm, ...
    'PVal', cfg.PVal, ...
    'CorrectionType', cfg.CorrectionType);

%% Plot (ITPC on top of power, same style as plotMultTransITPCEye)
fH = ITPCSimTest_plot(S, meta, 'CLim', cfg.CLim);

%% Optional: quick single-trial overlay figure
figure('Position',[880,50,500,350],'Color','w');
plot(t, sig(:,1:min(5,size(sig,2))));
hold on;
xline(0,'--k','LineWidth',1.5);
xlabel('sec (0 = absolute second 3)');
ylabel('uV');
title('Example trials (first up to 5)');
grid on;

fprintf('Done. Figure handle: fH\n');
