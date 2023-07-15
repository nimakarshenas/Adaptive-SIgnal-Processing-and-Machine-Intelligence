%% Part b - Implementing ALE
clc; clear all; close all;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

M = 5:5:20; 
delta = 1:25; step_size = 0.01; N=1000;
eta = filter([1, 0, 0.5], 1, randn(N+500, 1));
eta = eta(501:end);
x = sin((1:N)*0.01*pi)';
s = x + eta;

MSPE = zeros(length(M), length(delta));
L = length(s); 

figure(1);
subplot(1, 2, 1); hold on; 
set(gca,'fontsize', 16); 
for i = 1:length(M)
    for j = 1:length(delta)
        [~,xhat,~] = ale_lms(s, step_size, delta(j), M(i));
        MSPE(i, j) = (1/L)*(x-xhat)' * (x-xhat);
    end
    plot(delta, 10*log10(MSPE(i,:)), 'LineWidth', 1.5);
end
legend('M=5', 'M=10', 'M=15', 'M=20'); xlim([1, 25]);
xlabel('Delay ($\Delta$)');
grid on; grid minor;
ylabel('MSPE (dB)'); title('\textbf{MSPE vs. Filter Order \& Delay - ALE}');

[w, xhat, error] = ale_lms(s, 0.01, 5, 5);
subplot(1, 2, 2); hold on;
set(gca,'fontsize', 16); 
plot(1:length(s), s, 'LineWidth', 1.5);
plot(1:length(xhat), xhat, 'LineWidth', 1.5);
plot(1:length(x), x, 'LineWidth', 2); 
grid on; grid minor;
legend('Noisy signal', 'ALE; M=5, $\Delta=3$', 'Clean signal')
xlabel('Sample (n)'); ylabel('Signal'); title('\textbf{ALE Denoised Signal: Optimal Model}');
hold off;

%% Part b - model orders and delays
M = 5; delta = 3;
MSPE = zeros(2, 2);

figure(2);
% change model order fix delta
for i = 1:25
    [~,xhat,~] = ale_lms(s, step_size, delta, i);
    MSPE(1, i) = (1/L)*(x-xhat)' * (x-xhat);
    [~,xhat,~] = ale_lms(s, step_size, i, M);
    MSPE(2, i) = (1/L)*(x-xhat)' * (x-xhat);
end
subplot(1, 2, 1); hold on; set(gca,'fontsize', 16);
plot(1:25, 10*log10(MSPE(1,:)), LineWidth=1.5);  
hold off;
xlabel('Filter Order (M)'); ylabel('MSPE (dB)');
xlim([0 25]);
grid on; grid minor;
title('\textbf{MSPE vs. Filter Order} ($\Delta=3$)'); 

subplot(1, 2, 2); hold on; set(gca,'fontsize', 16); 
plot(1:25, 10*log10(MSPE(2,:)), LineWidth=1.5); 
hold off;
xlabel('Delay ($\Delta$)'); 
ylabel('MSPE (dB)');
grid on; grid minor;
xlim([0 25]);
title('\textbf{MSPE vs. Delay} (M=5)'); 

%% ale_lms function
function [w, xhat, error] = ale_lms(x, step_size, delta, M)

    w = zeros(M, length(x)+1); 
    error = zeros(size(x));
    xhat = x(size(x));
    
    for n = delta+M:length(x)
        u = flip(x(n-delta-M+1:n-delta));
        xhat(n) = dot(w(:, n), u);
        error(n) = x(n)  - xhat(n);
        w(:, n+1) = w(:, n) + step_size*error(n)*u;
    end
end