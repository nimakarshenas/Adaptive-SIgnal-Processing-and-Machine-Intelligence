clc; clear all; close all;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

sigma_n = 0.25;
data = filter([1], [1, -0.1, -0.8], sqrt(sigma_n)*randn(1500, 100));
data = data(501:end, :);

%% 1 realisation
[params1, error1] = lms_estimator(data(:, 1), 2, 0.05);
[params2, error2] = lms_estimator(data(:, 1), 2, 0.01);
figure(1); subplot(1, 2, 1); hold on; set(gca,'fontsize', 14);
plot(1:length(error1), 10*log10(error1.^2));
plot(1:length(error2), 10*log10(error2.^2));
xlabel('Sample'); ylabel('Squared Error (dB)');
title('\textbf{LMS Learning Curve (1 realisation)}');
legend('$\mu=0.05$','$\mu=0.01$');
ylim([-60 10])
grid on;
grid minor;
hold off;

%% 100 realisations
errors = zeros(2, 100, 1000); 
lrs = [0.05, 0.01];
for i = 1:100
    for j = 1:2
        [a, errors(j, i, :)] = lms_estimator(data(:, i), 2, lrs(j)); 
    end
end
errors = squeeze(mean(10*log10(errors.^2), 2));
figure(1); 
subplot(1, 2, 2); 
hold on; 
set(gca,'fontsize', 14);
plot(1:length(errors(1, :)), errors(1, :));
plot(1:length(errors(2, :)), errors(2, :));
xlabel('Sample'); ylabel('Squared Error (dB)');
title('\textbf{LMS Learning Curve (100 realisations)}');
legend('$\mathbf{\mu=0.05}$','$\mathbf{\mu=0.01}$', 'FontSize', 18);

grid on;
grid minor;
hold off;


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
