%% Part a)
clc; clear all; close all;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load PCAPCR.mat;
[U_x, S_x, V_x] = svd(X); 
[U_x_noise, S_x_noise, V_x_noise] = svd(Xnoise);

% Plotting singular values of X and Xnoise
figure(1); 
subplot(1, 2, 1); 
hold on;
stem(diag(S_x_noise), 'LineWidth', 1.5);
stem(diag(S_x), 'LineWidth', 1.5); 
grid on 
grid minor
set(gca,'fontsize', 13);
xlabel('Input variable index'); 
ylabel('Singular Value');
title('\textbf{Singular Values (SVs) of $\mathbf{X}$ and $\mathbf{X}_{noise}$ }'); 
xlim([1 10])
legend('$\mathbf{X}_{noise}$', '$\mathbf{X}$'); 
hold off;

% Square error between SVs
subplot(1, 2, 2);
stem((diag(S_x) - diag(S_x_noise)).^2, 'LineWidth', 1.5); 
set(gca,'fontsize', 13);
xlabel('Input variable index'); 
ylabel('Square Error');
title('\textbf{Square Error between SVs of $\mathbf{X}$ and $\mathbf{X}_{noise}$ }');
xlim([1 10])
grid on 
grid minor


%% Part b)
[U, S, V] = svds(Xnoise, 3); % keep the top 3 PCs
X_noise_lowrank = U*S*V';

% Plotting MSE
figure(2); 
hold on;
stem(mean(((X - Xnoise).^2), 1), 'LineWidth', 1.5); 
stem(mean(((X - X_noise_lowrank).^2), 1), 'LineWidth', 1.5); 
set(gca,'fontsize', 12); 
xlabel('Input variable index'); 
ylabel('MSE')
title('\textbf{MSE between $\mathbf{X}$ and \{ $\mathbf{X}_{noise}$, $\tilde{\mathbf{X}}_{noise}$}\} respectively');
legend('$\mathbf{X}_{noise}$', '$\tilde{\mathbf{X}}_{noise}$', 'Interpreter','latex', fontsize=14);
xlim([1 10])
grid on 
grid minor

hold off;


%% Part c)
% ols
B_ols = inv(Xnoise' * Xnoise)*X'*Y;

% pcr
[U, S, V] = svds(Xnoise, 3); 
for i = 1:length(diag(S)) 
    S(i, i) = 1/S(i, i);
end
B_pcr = (V * S' * U') * Y;


% Get output estimates
Y_ols = Xnoise * B_ols;
Y_pcr = Xnoise * B_pcr;

Ytest_ols = Xtest * B_ols;
Ytest_pcr = Xtest * B_pcr;

% squared training error
figure(4); 
subplot(1, 3, 1); 
hold on;
stem(sum((Y-Y_ols).^2, 1), 'LineWidth', 1.5); 
stem(sum((Y-Y_pcr).^2, 1), 'LineWidth', 1.5); 
set(gca,'fontsize', 14);
xlabel('Output Dimension'); 
ylabel('Squared Error');
title('\textbf{Squared Training Error}'); 
grid on;
grid minor;
legend('OLS', 'PCR'); hold off;

% Squared test error
subplot(1, 3, 2); hold on;
stem(sum((Ytest-Ytest_ols).^2, 1), 'LineWidth', 1.5); 
stem(sum((Ytest-Ytest_pcr).^2, 1), 'LineWidth', 1.5); 
set(gca,'fontsize', 14);
xlabel('Output Dimension'); 
ylabel('Squared Error');
grid on;
grid minor;
title('\textbf{Squared Test Error}'); 
legend('OLS', 'PCR'); 
hold off;


%% Part d)
ols_error = [];
pcr_error = [];

for i = 1:100
    [Y_ols, Y_ols_hat] = regval(B_ols);
    [Y_pcr, Y_pcr_hat] = regval(B_pcr);
    ols_error = [ols_error; sum((Y_ols - Y_ols_hat).^2, 1)];
    pcr_error = [pcr_error; sum((Y_pcr - Y_pcr_hat).^2, 1)];
end

ols_error = mean(ols_error, 1);
pcr_error = mean(pcr_error, 1);

subplot(1, 3, 3); 
hold on;
stem(ols_error, 'LineWidth', 1.5); 
stem(pcr_error, 'LineWidth', 1.5); 
set(gca,'fontsize', 14); 
xlabel('Output Dimension'); 
ylabel('Mean Squared Error')
title('\textbf{MSE of Estimate and Test Ensembles}');
legend('OLS', 'PCR'); 
grid on;
grid minor;
hold off;

