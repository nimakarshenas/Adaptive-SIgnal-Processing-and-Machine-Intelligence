clear all
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Load the sunspot data
load sunspot.dat
f_s=1; % Sampling frequency - One value per year

time = sunspot(:,1);
sunspot_series = sunspot(:,2);
L=length(sunspot_series);

%% Adding small a DC values to avoid zero values in log
sunspot_series=sunspot_series+eps;%distance from 1.0 to the next largest double-precision numbe

[psd_original_hann, f] = periodogram(sunspot_series,hann(L),L,f_s,'one-sided');
[psd_original_rect, f] = periodogram(sunspot_series,rectwin(L),L,f_s,'one-sided');


%% Removing mean and trend of series

%Remove mean only
sunspot_zero_mean = sunspot_series - mean(sunspot_series);

%Removing linear trend as well
sunspot_detrend = detrend(sunspot_zero_mean,1);

fft_detrend = fftshift(fft(sunspot_detrend,L));

[psd_detrend_hann, f] = periodogram(sunspot_detrend,hann(L),L,f_s,'one-sided');
[psd_detrend_rect, f] = periodogram(sunspot_detrend,rectwin(L),L,f_s,'one-sided');


%% Considering logarithm of data and removing mean
sunspot_log_zero_mean = log(sunspot_series)-mean(log(sunspot_series));

%Considering modified periodogram method (periodogram-based technique)
[psd_log_hann, f] = periodogram(sunspot_log_zero_mean,hann(L),L,f_s,'one-sided');
[psd_log_rect, f] = periodogram(sunspot_log_zero_mean,rectwin(L),L,f_s,'one-sided');


%% Now plotting the calculated values
figure(1);

% raw data power spectrum
subplot(1,3,1)
hold on
plot(f, psd_original_rect/10^6,'Linewidth',2)
plot(f, psd_original_hann/10^6,'Linewidth',2)

ylabel('\textbf{PSD$\times 10^6$}','FontSize',13)
xlabel('\textbf{Cycles/Year}','FontSize',13)
title('\textbf{PSD estimates of raw Sunspot data}','FontSize',13)
legend('\textbf{Rectanglur Window}', '\textbf{Hann Window}','FontSize',13)
grid on
grid minor
hold off

% zero-mean detrended data
subplot(1,3,2)
hold on
plot(f, psd_detrend_rect/10^6,'Linewidth',2)
plot(f, psd_detrend_hann/10^6,'Linewidth',2)
ylabel('\textbf{PSD$\times 10^6$}','FontSize',13)
xlabel('\textbf{Cycles/Year}','FontSize',13)
title({'\textbf{PSD estimates of zero-mean}','\textbf{+ detrended Sunspot data}'},'FontSize',13)
legend('\textbf{Rectanglur Window}', '\textbf{Hann Window}','FontSize',13)
grid on
grid minor
hold off

% log and zero-meaned data
subplot(1,3,3)
hold on
plot(f, psd_log_rect,'Linewidth',2)
plot(f, psd_log_hann,'Linewidth',2)
ylabel('\textbf{PSD}','FontSize',13)
xlabel('\textbf{Cycles/Year}','FontSize',13)
title({'\textbf{PSD estimates of log}', '\textbf{+zero-mean Sunspot data}'},'FontSize',13)
legend('\textbf{Rectanglur Window}', '\textbf{Hann Window}','FontSize',13)
grid on
grid minor
hold off