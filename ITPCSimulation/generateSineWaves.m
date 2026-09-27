function [signals, frequencies, phases] = generateSineWaves(N_sin, f1, f2, p, p1, p2, t)
%GENERATESINEWAVES Generate sine waves with random frequencies and phases.
%
% Inputs:
%   N_sin - Number of sine waves
%   f1,f2 - Minimum and maximum frequencies (Hz)
%   p      - Base phase (radians)
%   p1,p2  - Minimum and maximum phase jitter (radians)
%   t      - Time vector (seconds)
%
% Outputs:
%   signals     - N_sin-by-length(t) matrix
%   frequencies - Random frequency assigned to each sine wave
%   phases      - Resulting phase (p + jitter) of each sine wave

    frequencies = f1 + (f2 - f1) .* rand(N_sin, 1);
    jitter      = p1 + (p2 - p1) .* rand(N_sin, 1);
    phases      = p + jitter;

    signals = sin(2*pi .* frequencies .* t(:).' + phases);
end