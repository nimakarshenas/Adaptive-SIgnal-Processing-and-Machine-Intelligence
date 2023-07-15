close all
clear all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% General parameters

N=2000;
M=100; 
f_s=50;
f = [1, 1.5, 2];
t = (0:N-1)/f_s;
freqs = (-(2*N-1)/2:(2*N-1)/2-1)*(f_s/(2*N-1));

signal=[0.8*sin(2*pi*t*f(1)) + 1*sin(2*pi*t*f(2)) + 1.2*sin(2*pi*t*f(3)) zeros(1, N-length(t))];


figure(1);
%% Begin plot of absolute psd of 100 realisations and their empirical mean 

subplot(2,2,1)
psd = zeros(M,2*N-1);
for i=1:M
    signal_with_noise = signal + randn(1,N);
    [acf,~] = xcorr(signal_with_noise, 'biased');
    psd(i,:) = real(fftshift(fft(ifftshift(acf))));
    plot(freqs,psd(i,:),'c','Linewidth',1);
    hold on
end
xlim([0.5 2.5])
avg_psd = mean(psd);
plot(freqs,avg_psd,'b','Linewidth',1);
xlabel('Frequency (Hz)','FontSize',11); ylabel('PSD','FontSize',11)
title({'\textbf{Realisations and Mean of PSD}'},'FontSize',11)
std_psd = std(psd);
grid on
grid minor

subplot(2,2,2)
plot(freqs,std_psd,'r','Linewidth',1);
xlabel('Frequency (Hz)','FontSize',11); 
ylabel('$\hat{\sigma}_{P(\omega)}$','FontSize',11)
title({'\textbf{Standard deviation of PSD}'},'FontSize',11)
xlim([0.5 2.5])
grid on
grid minor


subplot(2,2,3)
psd_dB = zeros(M,2*N-1);
for i=1:M
   psd_dB(i,:) = pow2db(psd(i,:));
   plot(freqs,psd_dB(i,:),'c','Linewidth',1); 
   hold on 
end
xlim([0.5 2.5])
avg_psd_dB = mean(psd_dB);
plot(freqs,avg_psd_dB,'b','Linewidth',1);
xlabel('Frequency (Hz)','FontSize',11); ylabel('PSD (dB)','FontSize',11)
title({'\textbf{Realisations and Mean of PSD in dB}'},'FontSize',11)
grid on
grid minor


subplot(2,2,4)
std_psd_dB = std(psd_dB);
plot(freqs,std_psd_dB,'r','Linewidth',1);
xlim([0.5 2.5])
xlabel('Frequency (Hz)','FontSize',11); ylabel('$\hat{\sigma}_{P(\omega)}$(dB)','FontSize',11)
title({'\textbf{Standard deviation of the PSD in dB}'},'FontSize',11)
grid on 
grid minor