function [signals, cleanSignals, noise, frequencies, phases] = ...
    generateNoisySineWaves(N_sin, f1, f2, p, p1, p2, SNR_dB, t)
%GENERATENOISYSINEWAVES Generate sine waves with additive Gaussian noise.
%
% Inputs:
%   N_sin  - Number of sine waves
%   f1,f2  - Minimum and maximum frequencies (Hz)
%   p       - Base phase (radians)
%   p1,p2   - Minimum and maximum phase jitter (radians)
%   SNR_dB  - Desired signal-to-noise ratio in dB
%   t       - Time vector (seconds)
%
% Outputs:
%   signals      - Final noisy signals [N_sin x length(t)]
%   cleanSignals - Clean sine waves
%   noise        - Added zero-mean Gaussian noise
%   frequencies  - Randomly selected frequencies
%   phases       - Resulting phases: p + jitter

    % Random frequency and phase for each sine wave
    frequencies = f1 + (f2 - f1) .* rand(N_sin, 1);
    jitter      = p1 + (p2 - p1) .* rand(N_sin, 1);
    phases      = p + jitter;

    % Generate clean sine waves
    cleanSignals = sin(2*pi .* frequencies .* t(:).' + phases);

    % Generate zero-mean, unit-power Gaussian noise
    noise = randn(size(cleanSignals));
    noise = noise - mean(noise, 2);
    noise = noise ./ sqrt(mean(noise.^2, 2));

    % Required noise power: P_noise = P_signal / 10^(SNR_dB/10)
    signalPower = mean(cleanSignals.^2, 2);
    noisePower  = signalPower ./ (10.^(SNR_dB/10));

    % Scale and add the noise
    noise   = noise .* sqrt(noisePower);
    signals = cleanSignals + noise;
end