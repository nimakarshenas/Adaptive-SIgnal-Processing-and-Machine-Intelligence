clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

f_0 = 50; 
f_s = 5000; 
V = 10; 
phi = pi/3;
N = 1000;
n = 1:N;
Va = [10 15 10 10];
Vb = [10 10 5 10];
Vc = [10 10 10 5];
delta_b = [0 0 0 0];
delta_c = [0 0 0 0];
v = zeros(length(Va), length(n));
rhos = zeros(length(Va), 1);
cols = ["#0072BD", 	"#D95319", 	"#EDB120",	"#7E2F8E"];

figure(1); subplot(1,2,1); 
for idx = 1:length(Va)
    A = (sqrt(6)/6)*(Va(idx) + Vb(idx)*exp(1i*delta_b(idx)) + Vc(idx)*exp(1i*delta_c(idx)));
    B = (sqrt(6)/6)*(Va(idx) + Vb(idx)*exp(-1i*(delta_b(idx)+(2*pi/3))) + Vc(idx)*exp(-1i*(delta_c(idx)-(2*pi/3))));
    v(idx, :) =  A*exp(1i*2*pi*(f_0/f_s)*n + phi) + ...
    B*exp(-1i*2*pi*(f_0/f_s)*n + phi);
    [rhos(idx), ~] = circularity(v(idx, :));
    plot(v(idx, :), 'Marker','o', 'MarkerEdgeColor', cols(idx), 'MarkerFaceColor', cols(idx), MarkerSize=4); 
    hold on;
end
set(gca,'fontsize', 14);
grid on; grid minor;
xlabel('$\mathcal{R}$'); ylabel('$\mathcal{I}$');
legend([strcat('$\rho$ = ',num2str(rhos(1))) ', Balanced'],[strcat('$\rho$ = ',num2str(rhos(2))) strcat(', $V_{a}$ = ',num2str(Va(2))) strcat(', $V_{b}$= ',num2str(Vb(2))) strcat(', $V_{c}$= ',num2str(Vc(2)))],[strcat('$\rho$ = ',num2str(rhos(3))) strcat(', $V_{a}$ = ',num2str(Va(3))) strcat(', $V_{b}$= ',num2str(Vb(3))) strcat(', $V_{c}$= ',num2str(Vc(3)))], [strcat('$\rho$ = ',num2str(rhos(4))) strcat(', $V_{a}$ = ',num2str(Va(4))) strcat(', $V_{b}$= ',num2str(Vb(4))) strcat(', $V_{c}$= ',num2str(Vc(4)))])
title('Complex System Voltages for $\Delta_b, \Delta_c = 0$', 'and changing $V_a$, $V_b$, and $V_c$');

Va = [10 10 10 10];
Vb = [10 10 10 10];
Vc = [10 10 10 10];
delta_b = [0 pi/3 0 pi/3];
delta_c = [0 0 pi/3 pi/3];
v = zeros(length(Va), length(n));
rhos = zeros(length(Va), 1);
cols = ["#0072BD", 	"#D95319", 	"#EDB120",	"#7E2F8E"];
subplot(1,2,2); 
for idx = 1:length(Va)
    A = (sqrt(6)/6)*(Va(idx) + Vb(idx)*exp(1i*delta_b(idx)) + Vc(idx)*exp(1i*delta_c(idx)));
    B = (sqrt(6)/6)*(Va(idx) + Vb(idx)*exp(-1i*(delta_b(idx)+(2*pi/3))) + Vc(idx)*exp(-1i*(delta_c(idx)-(2*pi/3))));
    v(idx, :) =  A*exp(1i*2*pi*(f_0/f_s)*n + phi) + ...
    B*exp(-1i*2*pi*(f_0/f_s)*n + phi);
    [rhos(idx), ~] = circularity(v(idx, :));
    plot(v(idx, :), 'Marker','o', 'MarkerEdgeColor', cols(idx), 'MarkerFaceColor', cols(idx), MarkerSize=4); 
    hold on;
end
set(gca,'fontsize', 14);
grid on; grid minor;
xlabel('$\mathcal{R}$'); ylabel('$\mathcal{I}$');
legend([strcat('$\rho$ = ',num2str(rhos(1))) ', Balanced'],[strcat('$\rho$ = ',num2str(rhos(2))) ', $\Delta_{b} = \pi/3$, $\Delta_{c} = 0$'],[strcat('$\rho$ = ',num2str(rhos(3))) ', $\Delta_{b} = 0$, $\Delta_{c} = \pi/3$'],[strcat('$\rho$ = ',num2str(rhos(4))) ', $\Delta_{b} = \pi/3$, $\Delta_{c} = \pi/3$'])
title('Complex System Voltages for $V_a$, $V_b$, and $V_c$ = 10 ', 'and changing$\Delta_b, \Delta_c$');

%% functions
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