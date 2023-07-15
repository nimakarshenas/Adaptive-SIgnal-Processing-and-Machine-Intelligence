clc; clear all; close all;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

a1 = 0.1; 
a2 = 0.8; 
var = 0.25;
% realisation length
N = 2000;
% no. realisations
M = 100;
data = filter(1, [1, -0.1, -0.8], sqrt(var)*randn(N, M));
% remove transients 
data = data(501:end, :); 

% repeat for two step sizes
all_params = zeros(2, 2, 100, N-500); 
step_sizes = [0.05, 0.01];

for i = 1:M
    for j = 1:2
        [all_params(j, :, i, :), ~] = lms_estimator(data(:, i), 2, step_sizes(j)); 
    end
end
% take ensemble parameter estimate average at each time step
all_params = squeeze(mean(all_params, 3));
params_mu_005 = squeeze(all_params(1,:,:,:));
params_mu_001 = squeeze(all_params(2,:,:,:));

figure(1); 
subplot(1, 2, 1); 
hold on; 
set(gca,'fontsize', 15); 
ylim([0, 1]); 
plot(1:N-500, params_mu_005(1, :), 'LineWidth', 2);
plot(1:N-500, params_mu_005(2, :), 'LineWidth', 2);
plot(1:N-500, a1*ones(1,length(params_mu_005(1, :))), '--', 'LineWidth', 2);
plot(1:N-500, a2*ones(1,length(params_mu_005(2, :))), '--', 'LineWidth', 2);
xlabel('Time Step'); 
ylabel('Magnitude');
title('\textbf{LMS AR Estimation ($\mu=0.05$)}');
legend('$\hat{a}_{1}$','$\hat{a}_{2}$','$a_{1}$', '$a_{2}$');
grid on;
grid minor;
hold off;

subplot(1, 2, 2); 
hold on; 
set(gca,'fontsize', 15); 
ylim([0, 1]);
plot(1:N-500, params_mu_001(1, :), 'LineWidth', 2);
plot(1:N-500, params_mu_001(2, :), 'LineWidth', 2);
plot(1:N-500, a1*ones(1,length(params_mu_001(1, :))), '--', 'LineWidth',2);
plot(1:N-500, a2*ones(1,length(params_mu_001(2, :))), '--', 'LineWidth', 2);
xlabel('Time Step'); 
ylabel('Magnitude');
title('\textbf{LMS AR Estimation ($\mu=0.01$)}');
legend('$\hat{a}_{1}$','$\hat{a}_{2}$', '$a_{1}$', '$a_{2}$');
grid on;
grid minor;
hold off;

a_2_hat_001 = params_mu_001(2, end);
a_1_hat_001 = params_mu_001(1, end);

a_2_hat_005 = params_mu_005(2, end);
a_1_hat_005 = params_mu_005(1, end);



%% function for LMS estimation
function [params, error] = lms_estimator(data, order, step_size)

    params = zeros(order, length(data));
    error = zeros(size(data)) ;

    for i = order+1:length(data)
        current_error = data(i) - dot(flip(data(i-order:i-1)), params(:, i-1));
        error(i) = current_error;
        params(:, i) = params(:, i-1) + step_size*(current_error)*flip(data(i-order:i-1));
    end
end