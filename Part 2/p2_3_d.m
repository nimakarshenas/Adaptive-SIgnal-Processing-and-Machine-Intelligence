%% part d) removing mains component.
clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
load ../EEG_Data/EEG_Data_Assignment2.mat;

data = detrend(Cz);

n = 0:length(data)-1;
step_sizes = [0.1, 0.01, 0.001];
M = [5, 10, 20];
var = 0.1;

% plot original signal's periodogram
figure(1);
window_length = fs;
set(gca, "FontSize", 15, "FontWeight", "bold")
spectrogram(data, rectwin(window_length), round(0.2*(window_length)), 5*window_length, fs, 'yaxis');
title("\textbf{Cz with mains interference: Spectrogram}", FontSize=15)
ylim([0, 55]);
 % Cz with mains interference
mains = sin((2*pi*50/fs)*n) + sqrt(var)*randn(1, length(data));
mains = mains';
exp_no = 1;

figure(2);
for j = 1:length(step_sizes)
    for k = 1:length(M)
        [w,xhat] = anc_lms(data, mains, step_sizes(j), M(k));
        subplot(length(step_sizes), length(M), exp_no)
        ylim([0, 55]); 
        hold on; set(gca,'fontsize', 14);
        spectrogram(xhat, rectwin(window_length), round(0.2*(window_length)), 5*window_length, fs, 'yaxis');
        title(['M = ',num2str(M(k)),' and $\mu$ = ', num2str(step_sizes(j))], FontWeight="bold", FontSize=18);
        hold on;
        xlim([0 4.5])
        exp_no = exp_no + 1;
    end
end

[w,xhat] = anc_lms(data, mains, 0.01, 10);
[w,xhat2] = anc_lms(data, mains, 0.01, 20);
figure(3);
L = length(data);
[psd_x, f] = pwelch(data,hann(10000), 1000, 10000, fs,'one-sided');
[psd_xhat, f] = pwelch(xhat,hann(10000), 1000, 10000, fs,'one-sided');
[psd_xhat2, f] = pwelch(xhat2,hann(10000), 1000, 10000, fs,'one-sided');
plot(f, 10*log10(psd_x), LineWidth=1.5)
hold on;
plot(f, 10*log10(psd_xhat), LineWidth=1.5)
hold on;
plot(f, 10*log10(psd_xhat2), LineWidth=1.5)
grid on; grid minor;
set(gca, "FontSize", 16)
xlim([0 55])
title("\textbf{Periodogram of Cz and denoised Cz}")
legend("Original Signal", "ANC denoised: $\mu=0.01$, M=10", "ANC denoised: $\mu=0.01$, M=20")
ylabel("PSD (dB)")
xlabel("Frequency (Hz)")


%% anc lms function
function [w, xhat] = anc_lms(x, secondary_noise, step_size, M)
    
    w = zeros(M, length(x)); 
    eta = zeros(size(x));
    xhat = zeros(size(x));
    u = delayseq(repmat(secondary_noise, 1, M), [0:M-1])';
 
    for n = 1:length(x)
        eta(n) = dot(w(:, n), u(:, n));
        xhat(n) = x(n)  - eta(n);
        w(:, n+1) = w(:, n) + step_size*xhat(n)*u(:, n);
    end
end