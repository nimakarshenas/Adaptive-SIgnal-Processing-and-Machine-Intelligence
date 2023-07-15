%% 1.1. Properties of Power Spectral Density (PSD) 
clear all;
close all;
clc;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


N = 3000; 
f_s = 1000; 

t = 0:N-1;


%% calculations of P(\omega) related to unit impulse %%
delta_fn = zeros(1, N);
delta_fn(N/2) = 1;

% PSD from Equation 1
r_xx1 = xcorr(delta_fn, 'biased');

P_rxx1 = abs(fftshift(fft(r_xx1)));

% PSD from Equation 2
periodogram = (abs(fftshift(fft(delta_fn))).^2)./N;

% frequency values for Equation 1
f1 = -f_s/2:f_s/length(r_xx1):f_s/2 - f_s/length(r_xx1);
% frequency values for Equation 2
f2 = -f_s/2:f_s/length(periodogram):f_s/2 - f_s/length(periodogram);

%% calculations of P(\omega) related to Sinusoidally varying signal
sin_fn = sin(2*pi*0.5*t/f_s);
rxx_sin = xcorr(sin_fn, 'biased');

% PSD from Equation 1
P_rxx2 = abs(fftshift(fft(rxx_sin)));
% PSD from Equation 2
periodogram2 = abs(fftshift(fft(sin_fn))).^2 / N;


%% Plots related to unit impulse
figure(1);
subplot(3, 2, 1);
hold on;
plot(t, delta_fn);
set(gca,'fontsize', 16); 
xlabel('time (s)'); 
ylabel('Y(t)'); 
title('Unit Impulse');
grid on
grid minor
hold off;

subplot(3, 2, 3);
hold on;
plot((-N+1:N-1), r_xx1, 'r');
set(gca,'fontsize', 14); 
xlabel('Lag ($\tau$)'); 
ylabel('ACF'); 
title('ACF of unit impulse');
grid on
grid minor
hold off;

subplot(3, 2, 5); 
hold on;
plot(f1, P_rxx1, 'r'); 
plot(f2, periodogram, '--b');
set(gca,'fontsize', 14); 
xlabel('Frequency (Hz)'); 
ylabel('PSD'); 
title('PSD of Unit Impulse');
legend('Eq. (1)', 'Eq. (2)'); 
ylim([0 0.00033333333])
grid on
grid minor 
hold off;

%% Now plotting graphs related to sinusoid 

subplot(3, 2, 2);
hold on;
plot(t,sin_fn); 
set(gca,'fontsize', 14); 
xlabel('time (s)'); 
ylabel('Y(t)'); 
title('Sinusoidal signal');
grid on
grid minor
hold off;


subplot(3, 2, 4);
hold on;
plot((-N+1:N-1), rxx_sin);
set(gca,'fontsize', 14); 
xlabel('Lag ($\tau$)'); 
ylabel('ACF'); 
grid on
grid minor
title('ACF of Sinusoidal signal');
hold off;


subplot(3, 2, 6); 
hold on;
plot(f1, P_rxx2, 'r');
plot(f2, periodogram2, '--b');
set(gca,'fontsize', 14); 
xlabel('Frequency (Hz)'); 
ylabel('PSD'); 
title('PSD of Sinusoidal signal');
legend('Eq. (1)', 'Eq. (2)');
grid on
grid minor
hold off;