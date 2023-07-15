%% Part c - ANC vs ALE
clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

step_size = 0.01; N=1000;
eta = filter([1, 0, 0.5], 1, randn(N+500, 1));
eta = eta(501:end);
x = sin((1:N)*0.01*pi)';
s = x + eta;
secondary_noise = 0.7*eta + 0.1*delayseq(eta, 1) + 0.1*ones(length(eta), 1);
L = length(s);

MSPE = zeros(2, 2);

[~,xhat] = ale_lms(s, step_size, 3, 5);
[~,xhat_anc] = anc_lms(s, secondary_noise, step_size, 10);

figure(1); subplot(1, 2, 1);
hold on; set(gca,'fontsize', 18);
plot(1:length(xhat), xhat, LineWidth=1.5);
plot(1:length(xhat_anc), xhat_anc, LineWidth=1.5);
plot(1:length(x), x, 'LineWidth', 2);
xlabel('Sample (n)'); ylabel('Value'); 
title('ALE vs ANC: Signal Reconstruction of $x(n)$');
grid on; grid minor;
legend('ALE', 'ANC', '$x(n)$');
hold off;

MSPE_ale = 10*log10(movmean((x-xhat).^2, 40));
MSPE_anc = 10*log10(movmean((x-xhat_anc).^2, 40));

subplot(1, 2, 2);
hold on; set(gca,'fontsize', 18);
plot(1:length(MSPE_ale), MSPE_ale, LineWidth=1.5);
plot(1:length(MSPE_anc), MSPE_anc, LineWidth=1.5);
xlabel('Sample (n)'); ylabel('MSPE (dB)'); 
title('ALE vs ANC: MSPE');
legend('ALE', 'ANC');
grid on; grid minor;
hold off;

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