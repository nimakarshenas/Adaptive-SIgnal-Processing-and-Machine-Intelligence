close all
clear all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Generation of signals
L=10000;
t = (0:L-1);
f_s=1; %Sampling frequency for normalisation

%WGN
wgn = randn(L,1); %White Gaussian noise generated

%Noisy sinusoidal signal
f=0.4;
noisy_sine_wave = sin((2*pi*f).*t/f_s) + 4*randn(1,L);

%Filtered WGN - by MA filter
windowSize = 4; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
filtered_noise= filter(b,a,wgn);

%% Calculating the biased and unbiased ACF estimates

%WGN
[acf_WGN_biased, lag_WGN] = xcorr(wgn, 'biased'); %obtaining ACF estimate of WGN
[acf_WGN_unbiased,~] = xcorr(wgn, 'unbiased');


%Noisy sinusoidal
[acf_sine_biased, lag_sine] = xcorr(noisy_sine_wave, 'biased'); %obtaining ACF estimate of WGN
[acf_sine_unbiased,~] = xcorr(noisy_sine_wave, 'unbiased');

%Filtered WGN
[acf_filtered_biased, lag_filtered] = xcorr(filtered_noise, 'biased'); %obtaining ACF estimate of WGN
[acf_filtered_unbiased,~] = xcorr(filtered_noise, 'unbiased');


%% Plotting ACF estimates
figure(1);

subplot(3,2,1)
plot(lag_WGN,acf_WGN_unbiased,'Linewidth',1.5);
hold on
plot(lag_WGN,acf_WGN_biased,'Linewidth',1.5);
xlabel('\textbf{lag (k)}','FontSize',11)
ylabel('\textbf{ACF}','FontSize',11)
title('\textbf{ACF for WGN}','FontSize',11)
grid on 
grid minor



figure(1);
subplot(3,2,3)
plot(lag_sine,acf_sine_unbiased,'Linewidth',1.5);
hold on
plot(lag_sine,acf_sine_biased,'Linewidth',1.5);
xlabel('\textbf{lag (k)}','FontSize',11)
ylabel('\textbf{ACF}','FontSize',11)
legend('\textbf{Unbiased}','\textbf{Biased}','FontSize',10)
title('\textbf{ACF for noisy sinusoidal signal}','FontSize',11)
grid on 
grid minor

subplot(3,2,5)
plot(lag_filtered,acf_filtered_unbiased,'Linewidth',1.5);
hold on
plot(lag_filtered,acf_filtered_biased,'Linewidth',1.5);
xlabel('\textbf{lag (k)}','FontSize',11)
ylabel('\textbf{ACF}','FontSize',11)
title('\textbf{ACF for filtered WGN}','FontSize',11)
grid on 
grid minor


%% Obtaining the correlogram spectral estimator
% ifft and fft to bring ACF(0) to the start of the array containing ACF
% values, then taking the fftshift to obtain PSD

%For WGN
psd_wgn_unbiased = fftshift(fft(ifftshift(acf_WGN_unbiased)));
psd_wgn_biased =fftshift(fft(ifftshift(acf_WGN_biased)));

%For Noisy sinusoidal signal
psd_sine_unbiased = fftshift(fft(ifftshift(acf_sine_unbiased)));
psd_sine_biased =fftshift(fft(ifftshift(acf_sine_biased)));

%For filtered WGN
psd_filtered_unbiased = fftshift(fft(ifftshift(acf_filtered_unbiased)));
psd_filtered_biased =fftshift(fft(ifftshift(acf_filtered_biased)));

N=length(psd_wgn_unbiased);
%Definition of frequency axis
freqs = (-N/2:N/2-1)*(f_s/N);


%% Plotting Correlelograms
subplot(3,2,2)
plot(freqs,real(psd_wgn_unbiased),'Linewidth',1.5)
hold on
plot(freqs,real(psd_wgn_biased),'Linewidth',1.5)
ylabel('\textbf{PSD}','FontSize',11)
xlabel('\textbf{Normalised frequency} ($\pi$ rad/sample)','FontSize',11);
title('\textbf{Correlogram for filtered WGN}','FontSize',11)
grid on 
grid minor


subplot(3,2,4)
plot(freqs,real(psd_sine_unbiased),'Linewidth',1.5)
hold on
plot(freqs,real(psd_sine_biased),'Linewidth',1.5)
ylabel('\textbf{PSD}','FontSize',11)
xlabel('\textbf{Normalised frequency} ($\pi$ rad/sample)','FontSize',11);
legend('\textbf{Unbiased}','\textbf{Biased}','FontSize',10)
title('\textbf{Correlogram for Noisy sinusoidal signal}','FontSize',11)
grid on 
grid minor

subplot(3,2,6)
plot(freqs,real(psd_filtered_unbiased),'Linewidth',1.5)
hold on
plot(freqs,real(psd_filtered_biased),'Linewidth',1.5)
ylabel('\textbf{PSD}','FontSize',11)
xlabel('\textbf{Normalised frequency} ($\pi$ rad/sample)','FontSize',11);
title('\textbf{Correlogram for filtered WGN}','FontSize',11)
grid on 
grid minor