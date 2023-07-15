clc;
clear all;
close all;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

a1 = 0.1; 
a2 = 0.8; 
var = 0.25;
N=2000;
M=100;
data = filter(1, [1, -0.1, -0.8], sqrt(var)*randn(N, M));
data = data(501:end, :); 
step_sizes = [0.05, 0.01]; 
gammas = [0.1, 0.5, 0.9];
all_params = zeros(2, N-500, 100); 
figure(1);

%% mu = 0.05
for j = 1:3
    for i = 1:100
        [all_params(:, :, i), e] = leaky_lms_estimator(data(:, i), 2, step_sizes(1), gammas(j)); 
    end
    params = squeeze(mean(all_params, 3));
    subplot(2, 2, 1); 
    hold on;
    plot(1:length(params(1, :)), params(1,:),'LineWidth', 1.5); 
    
    subplot(2, 2, 2);
    hold on;
    plot(1:length(params(2, :)), params(2, :),'LineWidth',1.5);
end
subplot(2, 2, 1)
plot(1:length(params(1, :)), a1*ones(1,length(params(1, :))), '--','LineWidth', 1);
set(gca,'fontsize', 16); 
title("\textbf{LMS estimation of} $a_1$;  $\mu = 0.05$")
xlabel('Time Step'); 
ylabel('Estimate');
grid on;
grid minor;
hold off

subplot(2, 2, 2)
plot(1:length(params(2, :)), a2*ones(1,length(params(2, :))), '--','LineWidth', 1);
set(gca,'fontsize', 16); 
title("\textbf{LMS estimation of} $a_2$;  $\mu = 0.05$")
xlabel('Time Step'); 
ylabel('Estimate');
grid on;
grid minor;
hold off

%% mu = 0.01

for j = 1:3
    for i = 1:100
        [all_params(:, :, i), e] = leaky_lms_estimator(data(:, i), 2, step_sizes(2), gammas(j)); 
    end
    params = squeeze(mean(all_params, 3));
    subplot(2, 2, 3); 
    hold on;
    plot(1:length(params(1, :)), params(1,:),'LineWidth', 1.5); 
    
    subplot(2, 2, 4);
    hold on;
    plot(1:length(params(2, :)), params(2, :),'LineWidth',1.5);
end
subplot(2, 2, 3)
plot(1:length(params(1, :)), a1*ones(1,length(params(1, :))), '--','LineWidth', 1);
set(gca,'fontsize', 16); 
title("\textbf{LMS estimation of} $a_1$;  $\mu = 0.01$")
xlabel('Time Step'); 
ylabel('Estimate');
legend('$\gamma=0.1$','$\gamma=0.5$', '$\gamma=0.9$', 'True Value', 'Numcolumns', 2);
grid on;
grid minor;
hold off

subplot(2, 2, 4)
plot(1:length(params(2, :)), a2*ones(1,length(params(2, :))), '--','LineWidth', 1);
set(gca,'fontsize', 16); 
title("\textbf{LMS estimation of} $a_2$;  $\mu = 0.01$")
xlabel('Time Step'); 
ylabel('Estimate');
grid on;
grid minor;
hold off


%% function for leaky LMS estimation
function [params, error] = leaky_lms_estimator(data, order, step_size, gamma)

params = zeros(order, length(data)); 
error = zeros(size(data)) ;

for i = order+1:length(data)
    current_error = data(i) - dot(flip(data(i-order:i-1)), params(:, i-1));
    error(i) = current_error;
    params(:, i) = (1-step_size*gamma)*params(:, i-1) + step_size*(current_error)*flip(data(i-order:i-1));
end
end