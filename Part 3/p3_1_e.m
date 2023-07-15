clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

f_0 = 50; 
f_s = 1000; 
V = 10; 
phi = pi/3;
N = 200;
n = 1:N;
Va = [10 5 10];
Vb = [10 10 10];
Vc = [10 10 10];
delta_b = [0 0 pi/3];
delta_c = [0 0 0];
step_sizes_clms = [0.0001, 0.0001, 0.0001];
step_sizes_aclms = [0.0001, 0.0001, 0.0001];
v = zeros(length(Va), length(n));
rhos = zeros(length(Va), 1);
cols = ["#0072BD", 	"#D95319", 	"#EDB120",	"#7E2F8E"];
figure(1);
for idx = 1:length(Va)
    A = (sqrt(6)/6)*(Va(idx) + Vb(idx)*exp(1i*delta_b(idx)) + Vc(idx)*exp(1i*delta_c(idx)));
    B = (sqrt(6)/6)*(Va(idx) + Vb(idx)*exp(-1i*(delta_b(idx)+(2*pi/3))) + Vc(idx)*exp(-1i*(delta_c(idx)-(2*pi/3))));
    v(idx, :) =  A*exp(1i*2*pi*(f_0/f_s)*n + phi) + ...
    B*exp(-1i*2*pi*(f_0/f_s)*n + phi);
    [rhos(idx), ~] = circularity(v(idx, :));
    input_voltage = delayseq(v(idx, :)',1); 
    [h_clms, err_clms, ~] = clms(v(idx, :)',input_voltage,0,step_sizes_clms(idx),0);
    [q_aclms, err_aclms, ~] = aclms(v(idx,:)',input_voltage,0,step_sizes_aclms(idx),0);
    f_aclms = hg_to_f0(q_aclms,f_s);
    f_clms = h_to_f0(h_clms,f_s);
    
    subplot(1, 3, idx);
    plot([1:length(f_clms)-2], f_clms(3:end),'LineWidth',1.5);
    hold on
    plot([1:length(f_aclms)-2],f_aclms(3:end),'LineWidth',1.5);
    hold off
    yticks([0:25:100]);

end

subplot(1,3,1); hold on; set(gca,'fontsize', 16);
yticks([0:25:100]);
xlabel('Time (samples)'); ylabel('Frequency (Hz)');
title('Balaced System','Interpreter','Latex'); 
grid on; grid minor;
legend('CLMS','ACLMS'); ylim([0,100]);

subplot(1,3,2); hold on; set(gca,'fontsize', 16);
yticks([0,25,50,75,100]);
xlabel('Time (samples)'); ylabel('Frequency (Hz)');
title('Unbalanced System: $V_a \neq V_b, V_c$','Interpreter','Latex');
grid on; grid minor;
legend('CLMS','ACLMS'); ylim([0,100]);

subplot(1,3,3); hold on; set(gca,'fontsize', 16);
yticks([0,25,50,75,100]);
xlabel('Time (samples)'); ylabel('Frequency (Hz)');
title('Unbalanced System: $\Delta_b \neq 0$','Interpreter','Latex'); 
grid on; grid minor;
legend('CLMS','ACLMS'); ylim([0,100]);


%% functions
function [params, error, y_hat] = aclms(y, x, model_order, step_size, leak)

    params = zeros(2*(model_order+1), length(x),'like',1i);  
    params(model_order+1, :) = 1;
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
    params = ones(model_order+1, length(x),'like',1i); 
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

function f = h_to_f0(h,fs)
    f = (fs/(2*pi))*atan2(imag(h),real(h));
end

function f = hg_to_f0(q,fs)
    h = q(1,:); 
    g = q(2,:);
    f = (fs/(2*pi))*...
        atan2(real((sqrt(imag(h).^2 - abs(g).^2))),real(h));
end
