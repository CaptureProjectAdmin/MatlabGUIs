function [sig, t, meta] = ITPCSimTest_generateSignals(varargin)
%ITPCSimTest_generateSignals  Theta bursts with phase jitter in a 6 s window.
%
% Signals are 6 s epochs with a theta (4-8 Hz) sine burst centered near
% second 3. Across trials, phase is jittered within a configurable range.
% Time is returned centered on the burst (t=0 at second 3) to match the
% event-locked convention in plotMultTransITPCEye.

p = inputParser;
addParameter(p,'Fs',250);                 % Neuropace-like sample rate
addParameter(p,'DurationSec',6);
addParameter(p,'OnsetSec',3);             % burst center in absolute time
addParameter(p,'NTrials',50);
addParameter(p,'FreqHz',6);               % carrier in theta band (scalar or [lo hi])
addParameter(p,'PhaseJitterRad',pi/3);    % uniform phase jitter +/- this value
addParameter(p,'BurstHalfWidthSec',0.75); % Hann half-width around onset
addParameter(p,'Amplitude',50);           % uV-ish amplitude
addParameter(p,'NoiseStd',10);            % additive Gaussian noise (uV)
addParameter(p,'Seed',1);
parse(p,varargin{:});
opt = p.Results;

rng(opt.Seed,'twister');

n = round(opt.DurationSec * opt.Fs);
t_abs = (0:n-1)' ./ opt.Fs;               % 0 .. DurationSec
t = t_abs - opt.OnsetSec;                 % centered: -Onset .. +(Duration-Onset)

% Per-trial frequency in theta (4-8 Hz)
if numel(opt.FreqHz) == 1
    f_trial = repmat(opt.FreqHz,1,opt.NTrials);
else
    f_trial = opt.FreqHz(1) + diff(opt.FreqHz) * rand(1,opt.NTrials);
end
f_trial = min(max(f_trial,4),8);          % clamp to theta

% Per-trial phase jitter
phi = (2*rand(1,opt.NTrials) - 1) * opt.PhaseJitterRad;

% Burst envelope (Hann) centered at OnsetSec
env = zeros(n,1);
half = opt.BurstHalfWidthSec;
inburst = abs(t) <= half;
env(inburst) = 0.5 * (1 + cos(pi * t(inburst) / half));

sig = zeros(n,opt.NTrials);
for k = 1:opt.NTrials
    carrier = sin(2*pi*f_trial(k)*t_abs + phi(k));
    sig(:,k) = opt.Amplitude * env .* carrier + opt.NoiseStd * randn(n,1);
end

meta = struct();
meta.Fs = opt.Fs;
meta.DurationSec = opt.DurationSec;
meta.OnsetSec = opt.OnsetSec;
meta.NTrials = opt.NTrials;
meta.FreqHz = f_trial;
meta.PhaseRad = phi;
meta.PhaseJitterRad = opt.PhaseJitterRad;
meta.BurstHalfWidthSec = opt.BurstHalfWidthSec;
meta.Amplitude = opt.Amplitude;
meta.NoiseStd = opt.NoiseStd;
meta.Seed = opt.Seed;
meta.t_abs = t_abs;
end
