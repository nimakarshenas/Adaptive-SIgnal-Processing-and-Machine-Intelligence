close all; 
clear all;
clc;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load EEG_Data\EEG_Data_Assignment1.mat


%trim to get one experiment
%%POz = ecg2(269902:364206, 1);

N=length(POz);
f_s = 1200;
T_s = 1/f_s;
POz = POz - mean(POz);

% We want 5 dft samples per Hz, as recommended
nfft=5*f_s;
% Get periodogram using Hann window
[psd_POz f] = periodogram(POz,hann(N),nfft,f_s,'onesided');
psd_POz = pow2db(psd_POz);

figure(2);
subplot(2,2,1)
plot(f,psd_POz,'Linewidth',1)
xlim([0 60])
ylim([-150 -80])
set(gca,'fontsize', 14);
xlabel('Frequency(Hz)')
ylabel('PSD (dB)')
grid on
grid minor
title({'\textbf{Periodgram of EEG:}','\textbf{Standard method}'})

%% Calculate and plot psd with 10s windowing
%Window size 10s.
n_10 = 10/T_s; %Convert 10 seconds into sample size

[psd_10s,f_10s] = pwelch(POz,hann(n_10),0,nfft*(ceil(n_10/N)),f_s,'onesided');
psd_10s = pow2db(psd_10s); %Convert to dB

subplot(2,2,2)
plot(f_10s,psd_10s,'linewidth',1)
xlim([0 60])
set(gca,'fontsize', 14);
xlabel('Frequency(Hz)')
ylabel('PSD (dB)')
grid on
grid minor
title({'\textbf{Periodgram of EEG:}','\textbf{10s Windowing}'})

%% calculate and plot psd with 5s windowing
n_5 = 5/T_s;%Convert 5 seconds into sample size
[psd_5s,f_5s] = pwelch(POz,hann(n_5),0,nfft*(ceil(n_5/N)),f_s,'onesided');
psd_5s = pow2db(psd_5s); %Convert to dB
%Plot
subplot(2, 2, 3)
plot(f_5s,psd_5s,'linewidth',1)
xlim([0 60])
set(gca,'fontsize', 14);
xlabel('Frequency(Hz)')
ylabel('PSD (dB)')
grid on
grid minor
title({'\textbf{Periodgram of EEG:}','\textbf{5s Windowing}'})


%% calculate and plot psd with 1s windowing
n_1 = 1/T_s;%Convert 5 seconds into sample size
[psd_1s,f_1s] = pwelch(POz,hann(n_1),0,nfft*(ceil(n_1/N)),f_s,'onesided');
psd_1s = pow2db(psd_1s); %Convert to dB
%Plot
subplot(2, 2, 4)
plot(f_1s,psd_1s,'linewidth',1)
xlim([0 60])
set(gca,'fontsize', 14);
xlim([0 60])
xlabel('Frequency(Hz)')
ylabel('PSD (dB)')
grid on
grid minor
title({'\textbf{Periodgram of EEG:}','\textbf{1s Windowing}'})

%% Standard vs 10s averaged
figure(3)

plot(f,psd_POz,'Linewidth',1)
hold on
plot(f_10s,psd_10s,'r','linewidth',1)
set(gca,'fontsize', 14);
xlim([0 60])
ylim([-150 -80])
xlabel('Frequency(Hz)')
ylabel('PSD (dB)')
legend('Standard PSD','$\Delta_t$ = 10s, averaged')
sgtitle('\textbf{Comparison of Standard PSD and, Windowed and Averaged methods}')
grid on
grid minor

