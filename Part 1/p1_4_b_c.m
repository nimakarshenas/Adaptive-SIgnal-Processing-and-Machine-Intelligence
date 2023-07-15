clear all; 
close all; 
clc;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

f_s = 1;
N = 1000;
a = [2.76 -3.81 2.65 -0.92];
var = 1;
model_orders = 2: 14;
M = length(model_orders);

ar_model = arima('AR', a, 'Variance', var, 'Constant', 0);
ar_signal = simulate(ar_model, N);
ar_signal = ar_signal(N/2 + 1: end);
N = length(ar_signal);

[h, f] = freqz(1, [1 -a], N, f_s);
psd_actual = abs(h) .^ 2;
var_estimate = zeros(M, 1);

psd_a_m = cell(M, 1);

for m = 1: M
    [a_estimate, var_estimate(m)] = aryule(ar_signal, model_orders(m));
    freq_resp_ar_m = freqz(sqrt(var_estimate(m)), a_estimate, N);
    psd_a_m{m} = abs(freq_resp_ar_m) .^ 2;
end

figure;
% PSD: ground truth vs AR estimation
subplot(2, 2, 1);
plot(f, pow2db(psd_actual), 'LineWidth', 2);
hold on;
plot(f, pow2db(psd_a_m{1}), 'LineWidth', 2);
hold on;
plot(f, pow2db(psd_a_m{3}), 'LineWidth', 2);
hold on;
plot(f, pow2db(psd_a_m{7}), 'LineWidth', 2);
grid on;
grid minor;
legend('Truth', 'm = 2', 'm = 4', 'm = 8', 'FontSize', 11);
title({'\textbf{PSD estimate of an AR(4) model by an AR(m) model}', '\textbf{signal length = 500}'});
xlabel('Normalised frequency ($\pi$ rad/sample)', 'FontSize', 11);
ylabel('PSD (dB)', 'FontSize', 11);

% variance (noise power)
subplot(2, 2, 2);
plot(model_orders, pow2db(var_estimate), 'LineWidth', 2);
grid on;
grid minor;
title('\textbf{Variance of prediction error against AR order, N=500}', 'FontSize', 11);
xlabel('AR order (m)', 'FontSize', 11);
ylabel('Variance (dB)', 'FontSize', 11);



N = 10000;

ar_signal = simulate(ar_model, N);
ar_signal = ar_signal(1: N - 500);
N = length(ar_signal);
[h, f] = freqz(1, [1 -a], N, f_s);
psd_actual = abs(h) .^ 2;

var_estimate = zeros(M, 1);
psd_a_m = cell(M, 1);
for m = 1: M
    [a_estimate, var_estimate(m)] = aryule(ar_signal, model_orders(m));
    freq_resp_estimate = freqz(sqrt(var_estimate(m)), a_estimate, N);
    psd_a_m{m} = abs(freq_resp_estimate) .^ 2;
end


% PSD: ground truth vs AR estimation
subplot(2, 2, 3);
plot(f, pow2db(psd_actual), 'LineWidth', 2);
hold on;
plot(f, pow2db(psd_a_m{1}), 'LineWidth', 2);
hold on;
plot(f, pow2db(psd_a_m{3}), 'LineWidth', 2);
hold on;
plot(f, pow2db(psd_a_m{8}), 'LineWidth', 2);
grid on;
grid minor;
legend('Truth', 'm = 2', 'm = 4', 'm = 8', 'FontSize', 11);
title({'\textbf{PSD estimate of an AR(4) model by an AR(m) model}', '\textbf{signal length = 9500}'}, 'FontSize', 11);
xlabel('Normalised frequency ($\pi$ rad/sample)', 'FontSize', 11);
ylabel('PSD (dB)', 'FontSize', 11);


% variance (noise power)
subplot(2, 2, 4);
plot(model_orders, pow2db(var_estimate), 'LineWidth', 2);
grid on;
grid on;
grid minor;
title('\textbf{Variance of prediction error against AR order, N=9500}', 'FontSize', 11);
xlabel('AR order', 'FontSize', 11);
ylabel('Variance (dB)', 'FontSize', 11);