clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
load('time-series.mat');

y = y - mean(y); 
order = 4; N=length(y);
x = [zeros(order,1); y]; 
[pr,error,y_hat] = lms_standard(y, x, order, 1e-5, 0);

figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y, LineWidth=1);
plot(1:length(y_hat), y_hat, LineWidth=1);
xlabel('Samples (n)');
ylabel('Y'); title('\textbf{AR(4) Time Series Prediction}');
legend('True Signal', 'AR(4)'); hold off;
grid on; grid minor;

subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:N, y, LineWidth=1.5);
plot(1:N, y_hat, LineWidth=1.5);
grid on; grid minor;
xlim([800,1000]);
xlabel('Samples (n)');
ylabel('Y'); title('\textbf{AR(4) Time Series Prediction (Segment)}');
legend('True Signal', 'AR(4)'); hold off;

MSE = (error(800:end)'*error(800:end))/200;
Rp = pow2db(var(y_hat(800:end))/var(error(800:end)));


%% LMS function
function [params, err, y_hat] = lms_standard(y, x, order, mu, gamma)
    N = length(y);
    err = zeros(N,1);
    y_hat = zeros(N,1);
    params = zeros(order+1,N);

    for n = 1:N
        y_hat(n) = params(:,n).'*[1;x(n+order-1:-1:n)];
        err(n) = y(n) - y_hat(n);
        if n < N
            params(:, n+1) = (1-mu*gamma)*params(:, n)+ mu*err(n)*[1;x(n+order-1:-1:n)];
        end
    end
end