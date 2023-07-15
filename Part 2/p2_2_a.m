%% GASS
clc; clear all; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Setup parameters 
a_1 = 0.9; 
var = 0.5;
input = sqrt(var)*randn(1500, 1000);
output = filter([1, a_1], [1], input);
output = output(501:end, :);
input = input(501:end, :);
rho = 0.001; alpha = 0.8; lr0 = [0.2, 0.1, 0];
gass = {'standard', 'standard', 'ben', 'af', 'mx'};

figure(1);
for k = 1:3
    params_tot = zeros(5, 1000, 1000);
    lrs = [0.01, 0.1, lr0(k), lr0(k), lr0(k)]; 
    for j = 1:5
        for i = 1:1000
            [a, params_tot(j, :, i), e] = ...
                lms_arma(output(:, i), input(:, i), 0, 1, lrs(j), gass{j}, rho, alpha, 0);
        end
        param_error = -(squeeze(mean(params_tot(j, :, :), 3) - a_1));
        subplot(1, 3, k);
        hold on; set(gca,'fontsize', 15);
        plot(1:length(param_error), 10*log10(param_error), LineWidth=1.5); 
        hold off
    end
    subplot(1, 3, k);
    title(sprintf('$\\mu_{t=0}=%.1f$',lr0(k)), 'FontWeight','bold');
    legend('$\mu=0.01$','$\mu=0.1$', 'Benveniste', 'Ang \& Farhang', 'Matthews \& Xie', 'Interpreter', 'Latex'); 
    ylabel('Weight Error');
    grid on;
    grid minor;
    xlabel('Time Step');
end
sgtitle("\textbf{GASS Weight Error Curves}")

function [ar_params, ma_params, error] = lms_arma(output, input, p, q, step_size, algorithm, rho, alpha, leak)

    params = zeros(p+q, length(output)); 
    phi = zeros(p+q, length(output));
    step_sizes = step_size*ones(size(output));
    error = ones(size(output));
    
    for i = max([p,q])+1:length(output)-1
        flipped_data = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        error(i) = output(i) - dot(flipped_data, params(:, i));
        params(:, i+1) = (1-leak*step_size)*params(:, i) + step_sizes(i)*(error(i))*flipped_data;
        if strcmp(algorithm, 'af')
            if i > max([p,q])+1
                phi(:, i) = alpha*phi(:, i-1) + error(i-1)*prev_aug_dat;
            end
            step_sizes(i+1) = step_sizes(i) + rho*error(i)*dot(flipped_data, phi(:, i));
        elseif strcmp(algorithm, 'ben')
            if i > max([p,q])+1
                phi(:, i) = (eye(length(prev_aug_dat))-step_sizes(i-1)*(prev_aug_dat(:)*prev_aug_dat(:).'))...
                    *phi(:, i-1) + error(i-1)*prev_aug_dat;
            end
            step_sizes(i+1) = step_sizes(i) + rho*error(i)*dot(flipped_data, phi(:, i));
        elseif strcmp(algorithm, 'mx')
            if i > max([p,q])+1
                phi(:, i) = error(i-1)*prev_aug_dat;
            end
            step_sizes(i+1) = step_sizes(i) + rho*error(i)*dot(flipped_data, phi(:, i));
        elseif strcmp(algorithm, 'standard')
        end
        prev_aug_dat = flipped_data;
    end
    ar_params = params(1:p, :); ma_params = params(p+1:q, :);
end