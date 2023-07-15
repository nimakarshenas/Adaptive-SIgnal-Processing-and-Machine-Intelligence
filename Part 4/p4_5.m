clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
load('time-series.mat');

order = 4;
N = 1000;
x = [zeros(order,1); y];
[pr,error,y_hat] = lms_aug(y, x, order, 1e-7, 3, 0, 1);
figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y, LineWidth=1.5);
plot(1:length(y_hat), y_hat, LineWidth=1);
xlabel('Sample (n)');
ylabel('Value'); title('\textbf{Scaled LMS-Tanh Prediction}');
ylim([-50 50]);
grid on; grid minor;
legend('True Signal', 'LMS + tanh'); hold on;
subplot(1,2,2); 
plot(pr', LineWidth=1.5)
grid on; grid minor;set(gca,'fontsize', 16);
xlabel('Sample (n)');
ylabel('Parameter Value'); title('\textbf{Evolution of Parameter values}');
legend('$b$','$w_1$','$w_2$','$w_3$','$w_4$'); hold on;


MSE = (error(800:end)'*error(800:end))/200;
Rp = pow2db(var(y_hat(800:end))/var(error(800:end)));
MSE_full = (error'*error)/N;
Rp_full = pow2db(var(y_hat)/var(error));


function [params, err, y_hat, a] = lms_aug(y, x, order, mu_w, mu_a, gamma, pt)
    N = length(y);
    err = zeros(N,1);
    y_hat = zeros(N,1);
    if pt
        [params, a] = pretrain(y,x,100,20,order, mu_w, mu_a, gamma);
    else
         params = zeros(order+1,N); 
         a = 80*ones(N,1);
    end

    for n = 1:N
        y_hat(n) = a(n)*params(:,n).'*[1; x(n+order-1:-1:n)];
        err(n) = y(n) - y_hat(n);
        if n < N
            params(:, n+1) = (1-mu_w*gamma)*params(:, n)+ mu_w*a(n)*err(n)*(1-(tanh(params(:,n).'*[1; x(n+order-1:-1:n)]).^2))*[1; x(n+order-1:-1:n)];
            a(n+1) = a(n) + mu_a*err(n)*tanh(params(:,n).'*[1; x(n+order-1:-1:n)]);
        end
    end
end

function [params, a] = pretrain(y,x,n_epochs,K,order, mu_w, mu_a, gamma)
    N = length(y);
    params = zeros(order+1,N);
    weights = zeros(order+1,1);
    err = zeros(N,1);
    y_hat = zeros(N,1);
    a = 80*ones(N,1);
    for e = 1:n_epochs
        for k = 1:K
            for n = 1:N
                y_hat(k) = a(k)*weights.'*[1; x(k+order-1:-1:k)];
                err(k) = y(k) - y_hat(k);
                if n < N
                    weights = (1-mu_w*gamma)*weights+ mu_w*a(k)*err(k)*(1-(tanh(weights.'*[1; x(k+order-1:-1:k)]).^2))*[1; x(k+order-1:-1:k)];
                    a(k+1) = a(k) + mu_a*err(k)*tanh(weights.'*[1; x(k+order-1:-1:k)]);
                end
            end
        end
    end
    params(:,1:K) = repmat(weights, 1, K);
end