clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
load('time-series.mat');

order = 4;
N = 1000;
x = [zeros(order,1); y];
[pr,error,y_hat] = lms_bias(y, x, order, 1e-7, 3, 0);
MSE_adap = (error(800:end)'*error(800:end))/200;
Rp_adap = pow2db(var(y_hat(800:end))/var(error(800:end)));
figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y, LineWidth=1.5);
plot(1:length(y_hat), y_hat, LineWidth=1);
xlabel('Sample (n)');
ylabel('Value'); title('\textbf{Scaled LMS-Tanh Prediction}');
ylim([-50 50]);
grid on; grid minor;
legend('True Signal', 'LMS + tanh'); hold on;
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y, LineWidth=1.5);
plot(1:length(y_hat), y_hat, LineWidth=1);
xlim([800,1000]);
grid on; grid minor;
xlabel('Sample (n)');
ylabel('Value'); title('\textbf{Scaled LMS-Tanh Prediction (Zoomed)}');
legend('True Signal', 'LMS + tanh'); hold on;

figure(2);
plot(pr', LineWidth=1.5)
grid on; grid minor;set(gca,'fontsize', 16);
xlabel('Sample (n)');
ylabel('Parameter Value'); title('\textbf{Evolution of Parameter values}');
legend('$b$','$w_1$','$w_2$','$w_3$','$w_4$'); hold on;



MSE = (error(800:end)'*error(800:end))/200;
Rp = pow2db(var(y_hat(800:end))/var(error(800:end)));

MSE_full = (error'*error)/N;
Rp_full = pow2db(var(y_hat)/var(error));


function [params, err, y_hat, a] = lms_bias(y, x, order, mu_w, mu_a, gamma)
    N = length(y);
    err = zeros(N,1);
    y_hat = zeros(N,1);
    params = zeros(order+1,N);
    a = 80*ones(N,1);

    for n = 1:N
        y_hat(n) = a(n)*params(:,n).'*[1; x(n+order-1:-1:n)];
        err(n) = y(n) - y_hat(n);
        if n < N
            params(:, n+1) = (1-mu_w*gamma)*params(:, n)+ mu_w*a(n)*err(n)*(1-(tanh(params(:,n).'*[1; x(n+order-1:-1:n)]).^2))*[1; x(n+order-1:-1:n)];
            a(n+1) = a(n) + mu_a*err(n)*tanh(params(:,n).'*[1; x(n+order-1:-1:n)]);
        end
    end
end