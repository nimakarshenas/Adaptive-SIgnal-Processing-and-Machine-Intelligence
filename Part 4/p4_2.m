clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
load('time-series.mat');

y = y - mean(y); 
order = 4;
x = [zeros(order,1); y]; 
[pr,error,y_hat] = lms_percept(y, x, order, 2e-7, 0);

figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y, LineWidth=1);
plot(1:length(y_hat), y_hat, LineWidth=1);
xlabel('Sample (n)');
ylabel('Value'); title('\textbf{LMS Prediction with Tanh Activation (Segment)}');
legend('True Signal', 'LMS + tanh'); hold on; grid on; grid minor;

subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y, LineWidth=1.5);
plot(1:length(y_hat), y_hat, LineWidth=1.5);
xlim([800,1000]);
xlabel('Sample (n)');
ylabel('Value'); title('\textbf{LMS Prediction with Tanh Activation (Segment)}');
legend('True Signal', 'LMS + tanh'); hold on; grid on; grid minor;

MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));


%% LMS 
function [params, err, y_hat] = lms_percept(y, x, order, mu, gamma)
    N = length(y);
    err = zeros(N,1);
    y_hat = zeros(N,1);
    params = zeros(order,N);

    for n = 1:N
        y_hat(n) = params(:,n).'*[x(n+order-1:-1:n)];
        err(n) = y(n) - y_hat(n);
        if n < N
            params(:, n+1) = (1-mu*gamma)*params(:, n)+ mu*err(n)*(1-(tanh(params(:,n).'*[x(n+order-1:-1:n)]).^2))*[x(n+order-1:-1:n)];
        end
    end
end