clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
load('time-series.mat');

y = y - mean(y); 
order = 4;
errors = ones(1, 62);
N = length(y);
no_experiments = 100;
a_min = 38;
a_max = 100;
a = a_min:(a_max-a_min)/no_experiments:a_max;
x = [zeros(order,1); y]; 
erros = zeros(1, no_experiments);
for i = 1:no_experiments
    [pr,error,y_hat] = lms_scale(y, x, order, 1e-7, a(i), 0);
    errors(i) = pow2db(mean(error(1:end).^2));
end
figure(1); 
subplot(1, 2, 1)
plot(a(1:end-1), errors, LineWidth=1); set(gca,'fontsize', 16);
xlabel('$a$');
ylabel('MSE(dB)'); title('\textbf{Finding optimal value for} $a$');
grid on; grid minor;
subplot(1, 2, 2); 
[pr,error,y_hat] = lms_scale(y, x, order, 1e-7,90.7, 0);
MSE_fixed = (error(800:end)'*error(800:end))/200;
Rp_fixed = pow2db(var(y_hat(800:end))/var(error(800:end)));

plot(1:length(y), y, LineWidth=1.5);
hold on;
plot(1:length(y_hat), y_hat, LineWidth=1);
set(gca,'fontsize', 16);
xlim([0,1000]);
grid on; grid minor;
xlabel('Sample(n)');
ylabel('Y'); title('\textbf{Scaled LMS-Tanh Prediction(Zoomed)}: $a_{opt}=96.342$');
legend('True Signal', 'LMS + tanh');

%% Dynamic 
y = y - mean(y); 
[pr,error,y_hat, a] = lms_dynamic_scale(y, x, order, 1e-7, 3, 0);

figure(2); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(a, LineWidth=1.5);
xlabel('Sample (n)');
ylabel('$a(n)$'); title('\textbf{Evolution of }$a$');
grid on; grid minor; hold on;
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y, LineWidth=1.5);
plot(1:length(y_hat), y_hat, LineWidth=1);
xlim([0,1000]);
grid on; grid minor;
xlabel('Sample (n)');
ylabel('Value'); title('\textbf{Scaled LMS-Tanh Prediction (Zoomed)}');
legend('True Signal', 'LMS + tanh'); hold on;

MSE = (error(800:end)'*error(800:end))/200;
Rp = pow2db(var(y_hat(800:end))/var(error(800:end)));
MSE_full = (error'*error)/N;
Rp_full = pow2db(var(y_hat)/var(error));

%% LMS 
function [params, err, y_hat] = lms_scale(y, x, order, mu_w, a, gamma)
    N = length(y);
    err = zeros(N,1);
    y_hat = zeros(N,1);
    params = zeros(order+1,N);

    for n = 1:N
        y_hat(n) = a*params(:,n).'*[1;x(n+order-1:-1:n)];
        err(n) = y(n) - y_hat(n);
        if n < N
            params(:, n+1) = (1-mu_w*gamma)*params(:, n)+ mu_w*a*err(n)*(1-(tanh(params(:,n).'*[1;x(n+order-1:-1:n)]).^2))*[1;x(n+order-1:-1:n)];
        end
    end
end

function [params, err, y_hat, a] = lms_dynamic_scale(y, x, order, mu_w, mu_a, gamma)
    N = length(y);
    err = zeros(N,1);
    y_hat = zeros(N,1);
    params = zeros(order,N);
    a = 80*ones(N,1);

    for n = 1:N
        y_hat(n) = a(n)*params(:,n).'*[x(n+order-1:-1:n)];
        err(n) = y(n) - y_hat(n);
        if n < N
            params(:, n+1) = (1-mu_w*gamma)*params(:, n)+ mu_w*a(n)*err(n)*(1-(tanh(params(:,n).'*[x(n+order-1:-1:n)]).^2))*[x(n+order-1:-1:n)];
            a(n+1) = a(n) + mu_a*err(n)*tanh(params(:,n).'*[x(n+order-1:-1:n)]);
        end
    end
end