clc; clear all; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

b1 = 1.5 + 1i; 
b2 = 2.5 - 0.5i;

N=1000;
var = 1;
noise_power = pow2db(var);
x = wgn(N,1,noise_power,'complex');
y = zeros(N,1,'like',1i);
error_aclms = zeros(N,1);
error_clms = zeros(N,1);

% generate WLSA
for n = 1:1000
    % need at least one sample
    if n == 1
        y(n) = 0i;
    else
        y(n) = x(n) + b1*x(n-1) + b2*conj(x(n-1));
    end
end
[~,error_aclms(:),~] = aclms(y(:),x(:),1, 0.1, 0);
[~,error_clms(:),~] = clms(y(:),x(:),1, 0.1, 0);

mean_error_aclms = mean(abs(error_aclms).^2,2);
mean_error_clms = mean(abs(error_clms).^2,2);
% Plots
figure(1); subplot(1,2,1); hold on; set(gca,'fontsize', 15);
plot(y, 'ro'); 
plot(x, 'bo'); 
title('\textbf{Circularity Plot}'); 
xlabel('$\mathcal{R}$'); ylabel('$\mathcal{I}$'); 
xlim([-15 15]); ylim([-15 15]);
legend("WLMA(1)", "WGN")
grid on; grid minor;
hold off;

% error curves
subplot(1,2,2); hold on; 
set(gca,'fontsize', 15);
plot([1:1000], pow2db(mean_error_aclms), LineWidth=1.5);
plot([1:1000], pow2db(mean_error_clms), LineWidth=1.5);
ylabel('Squared Error (dB)'); xlabel('Time (samples)');

grid on; grid minor;
title('\textbf{Learning Curve}'); 
legend('ACMLS', 'CMLS');
hold off;

[eta_x, rho_x] = circularity(x);
[eta_y, rho_y] = circularity(y);

function [params, error, y_hat] = aclms(y, x, model_order, step_size, leak)

    params = zeros(2*(model_order+1), length(x),'like',1i);  
    error = zeros(size(x),'like',1i);
    y_hat = zeros(size(x),'like',1i);
    x_pad = [zeros(model_order,1,'like',1i); x];
    for n = 1:length(x)
        x_aug = [x_pad(n+model_order:-1:n); conj(x_pad(n+model_order:-1:n))];
        y_hat(n) = params(:,n)'*x_aug;
        error(n) = y(n) - y_hat(n);
        if n < length(x)
            params(:, n+1) = (1-step_size*leak)*params(:, n) + step_size*conj(error(n))*x_aug;
        end
    end
end

function [params, error, y_hat] = clms(y, x, model_order, step_size, leak)
    params = zeros(model_order+1, length(x),'like',1i); 
    error = zeros(size(x),'like',1i);
    y_hat = zeros(size(x),'like',1i);
    x_pad = [zeros(model_order,1); x];
    for n = 1:length(x)
        y_hat(n) = params(:,n)'*x_pad(n+model_order:-1:n); 
        error(n) = y(n) - y_hat(n);
        if n < length(x)
            params(:, n+1) = (1-step_size*leak)*params(:, n) + step_size*conj(error(n))*x_pad(n+model_order:-1:n);
        end
    end
end

function [eta,rho] = circularity(data)
    num = mean(data.*data); den = mean(abs(data).^2);
    rho = num/den;
    eta = abs(rho);
end