clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

n_samples = 1500;
f = [100*ones(500,1); ...
    100*ones(500,1) + ([501:1000]' - 500)/2; ...
    100*ones(500,1) + (([1001:1500]' -  1000)/25).^2];
phi = cumsum(f);
f_s = 1500;
variance = 0.05;
eta = wgn(n_samples,1,pow2db(variance),'complex');
y = exp(1i*2*pi*phi/f_s) + eta;

input = delayseq(y,1);
lrs = [1e-1, 1e-2,1e-3];
figure(1);
for idx = 1:3
    [a,~,~] = clms(y, input, 0, lrs(idx), 0);
    H = zeros(n_samples,n_samples);
    for n = 1:n_samples
        [h, w] = freqz(1 , [1; -conj(a(n))], n_samples,f_s); 
        H(:, n) = abs(h).^2; 
    end

    medianH = 50*median(median(H));
    H(H > medianH) = medianH;
    subplot(1,3,idx); hold on; set(gca,'fontsize', 16);
    mesh([1:n_samples], w, H);
    view(2);
    xlabel('Time (samples)');
    ylabel('Frequency (Hz)');
    ylim([0,500]);
    title(strcat('CLMS Spectral Estimate ($\mu=',num2str(lrs(idx)),'$)'), 'Interpreter', 'Latex');
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