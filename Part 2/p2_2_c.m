clc; clear all; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

a_1 = 0.9; var = 0.5;
N = 1000;
input = sqrt(var)*randn(N, 10000);
output = filter([1, a_1], [1], input);
output = output(501:end, :);
input = input(501:end, :);

rho_gngd = 0.005; 
rho_bene = 0.002;
beta_gngd = 1; 
step_size_init_bene = 0.1;


params_tot_g = zeros(N-500, 10000); 
error_tot_g = zeros(N-500, 10000);
params_tot_b = zeros(N-500, 10000); 
error_tot_b = zeros(N-500, 10000);


for i = 1:10000
    if mod(i, 100) == 0
    end
    [~, params_tot_g(:, i), error_tot_g(:, i)] = ...
        gngd(output(:, i), input(:, i), 0, 1, beta_gngd, rho_gngd, 0);
    [~, params_tot_b(:, i), error_tot_b(:, i)] = ...
        lms_arma(output(:, i), input(:, i), 0, 1, step_size_init_bene,'ben', rho_bene, 0, 0);
end
params_tot_g = -(squeeze(mean(params_tot_g, 2)) - a_1);
params_tot_b = -(squeeze(mean(params_tot_b, 2)) - a_1);
error_tot_g = squeeze(mean(mag2db(error_tot_g.^2), 2));
error_tot_b = squeeze(mean(mag2db(error_tot_b.^2), 2));

figure(1);
subplot(1,2,1) 
set(gca,'fontsize', 16); 
hold on;
xlim([0 400])
plot(1:length(params_tot_g), params_tot_g, LineWidth=1.5);
plot(1:length(params_tot_b), params_tot_b, LineWidth=1.5);
title('\textbf{Weight Error}')
sgtitle('\textbf{GNGD vs Benveniste}', fontsize=16)
legend('\textbf{GNGD}', '\textbf{Benveniste}');
ylabel('Weight Error'); xlabel('Time Step');
grid on; grid minor;
hold off;

subplot(1,2,2)
xlim([0 80]); 
set(gca,'fontsize', 16); 
hold on;
plot(1:length(params_tot_g), params_tot_g, LineWidth=1.5);
plot(1:length(params_tot_b), params_tot_b, LineWidth=1.5);
title('\textbf{Weight Error (Zoomed)}')
legend('\textbf{GNGD}', '\textbf{Benveniste}');
ylabel('Weight Error'); xlabel('Time Step');
grid on; grid minor;
hold off;




%% LMS and GASS function
function [ar_params, ma_params, prediction_error] = lms_arma(output, input, p, q, step_size, algorithm, rho, alpha, leak)

    params = zeros(p+q, length(output)); 
    phi = zeros(p+q, length(output));
    step_sizes = step_size*ones(size(output));
    prediction_error = ones(size(output));
    
    for i = max([p,q])+1:length(output)-1
        req_data = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        prediction_error(i) = output(i) - dot(req_data, params(:, i));
        params(:, i+1) = (1-leak*step_size)*params(:, i) + step_sizes(i)*(prediction_error(i))*req_data;
        if strcmp(algorithm, 'af')
            if i > max([p,q])+1
                phi(:, i) = alpha*phi(:, i-1) + prediction_error(i-1)*prev_aug_dat;
            end
            step_sizes(i+1) = step_sizes(i) + rho*prediction_error(i)*dot(req_data, phi(:, i));
        elseif strcmp(algorithm, 'ben')
            if i > max([p,q])+1
                phi(:, i) = (eye(length(prev_aug_dat))-step_sizes(i-1)*(prev_aug_dat(:)*prev_aug_dat(:).'))...
                    *phi(:, i-1) + prediction_error(i-1)*prev_aug_dat;
            end
            step_sizes(i+1) = step_sizes(i) + rho*prediction_error(i)*dot(req_data, phi(:, i));
        elseif strcmp(algorithm, 'mx')
            if i > max([p,q])+1
                phi(:, i) = prediction_error(i-1)*prev_aug_dat;
            end
            step_sizes(i+1) = step_sizes(i) + rho*prediction_error(i)*dot(req_data, phi(:, i));
        elseif strcmp(algorithm, 'standard')
        end
        prev_aug_dat = req_data;
    end
    ar_params = params(1:p, :); ma_params = params(p+1:q, :);
end


%% GNGD function 
function [ar_params, ma_params, prediction_error] = gngd(output, input, p, q, beta, rho, leak)

    params = zeros(p+q, length(output));
    prediction_error = ones(size(output));
    eps = ones(size(output))/beta;
    
    for i = max([p,q])+1:length(output)-1
        req_data = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        prediction_error(i) = output(i) - dot(req_data, params(:, i));
        step_size = beta/(eps(i) + dot(req_data, req_data));
        params(:, i+1) = (1-leak*beta)*params(:, i) + step_size*(prediction_error(i))*req_data;
        % epsilon update requires at least one sample of past data
        if i>max([p,q])+1
            eps(i+1) = eps(i) - (rho*beta*prediction_error(i)*prediction_error(i-1)*dot(old_dat, req_data))/( eps(i-1) + dot(old_dat, old_dat) ).^2;
        end
        old_dat = req_data;
        
    end
    ar_params = params(1:p, :); ma_params = params(p+1:q, :);
end