close all
clear


Fs = 2500;
t = 1/Fs:1/Fs:10;
p1 = pi/4;
p2 = pi/2;
p = pi/3;
f1 = 30; % Frequency of the first sine wave
f2 = 30; % Frequency of the second sine wave
N_sin = 1;
SNR_dB = 90;

% [Y, f, phase] = generateSineWaves(N_sin,f1, f2, p, p1, p2, t);

[YN, Y, noise, f, phases] = ...
    generateNoisySineWaves(N_sin, f1, f2, p, p1, p2, SNR_dB, t);


figure
plot(t, Y);
xlabel('Time (s)');
ylabel('Amplitude');


[f,coi,cfs] = morseSpecGram(Y',Fs,[1,80]);

CF = squeeze(mean(abs(cfs),3));

figure
imagesc(t,f,CF')
axis xy
xlabel("Time (s)")
ylabel("Frequency (Hz)")


figure
spectrogram(Y,'yaxis')

