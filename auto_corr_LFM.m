clc; clear all; close all;
fc=20000; %centre frequency 20khz
Tc= 0.05;
miue= 100000;%Hz/s;
BT=miue*(Tc)^2; % Time-Bandwidth Product
BW=BT/Tc;% Bandwidth
Fstart=17500; %17.5 kHz
Fe=22500; %22.5 kHz
Fs=6*Fe;
Ts=1/Fs;
n=0:10000;
nTs=0:Ts:Tc;
LFM=cos(2*pi*fc.*nTs + miue*pi.*(nTs).^2).* rect(nTs,nTs/Tc);
figure;plot(nTs*10^3,LFM); title('LFM chirp signal')
xlabel('time (ms)');ylabel('amplitude (v)');
figure;
spectrogram(LFM,256,250,[],6*22500,'yaxis')
title('Spectrogram of transmitted LFM')
LFM1= flip(LFM);
% Bonus Point
%%autocorrelation 
yy=my_conv(LFM,LFM1);
time=((1:length(yy))/Fs)*1e3;
figure; plot(time,abs(yy)); title('auto-correlation of LFM chirp');
xlabel('time (ms)');ylabel('amplitude (v)');
