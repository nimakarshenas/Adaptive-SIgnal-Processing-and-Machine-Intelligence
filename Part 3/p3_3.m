clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


n_samples = 1500;
f_mag = [100*ones(500,1); ...
    100*ones(500,1) + ([501:1000]' - 500)/2; ...
    100*ones(500,1) + (([1001:1500]' -  1000)/25).^2];
phi = cumsum(f_mag);
fs = 1000;
variance = 0.05;
eta = wgn(n_samples,1,pow2db(variance),'complex');
y = exp(1i*2*pi*phi/fs) + eta;
t = 1:n_samples;
k = 0:fs-1;
x = (1/fs)*exp(1i*2*pi*k'*t./fs);

leaks = [0, 0.1, 0.3];
figure(1);
for idx = 1:3
    [w1,~] = dft_clms(y.', x, leaks(idx));
    subplot(1,3,idx); hold on; set(gca,'fontsize', 16);
    mesh(t, k(1:floor(fs/2)), abs(w1(1:floor(fs/2),:)).^2);
    view(2);
    xlabel('Time (samples)');
    ylabel('Frequency (Hz)');
    ylim([0,500]);
    title(strcat('CLMS Fourier Estimate ($\gamma=',num2str(leaks(idx)),'$)'),'Interpreter','Latex');
    hold off;
end

%% Part d

load('..\EEG_Data\EEG_Data_Assignment2.mat');

a = 1;
t_range = a:a+1199;
POz = detrend(POz(t_range));
n_samples = length(t_range);
f_ax = 0:fs-1;
x = (1/fs)*exp(1i*2*pi*f_ax'*t_range./fs);
[w1,~] = dft_clms(POz, x, 0);
figure(2); 
hold on; set(gca,'fontsize', 16);
mesh(t_range, f_ax(1:floor(fs/2)), abs(w1(1:floor(fs/2),:)).^2);
view(2); ylim([0,60]);
xlabel('Time (samples)');
ylabel('Frequency (Hz)');
title('CLMS Spectrogram Estimate (POz): $\gamma=0$');



function [w, error] = dft_clms(y, x, leak)

    [nft, n_samples] = size(x);
    w = zeros(nft, n_samples,'like',1i); 
    error = zeros(n_samples,1,'like',1i);
    mu = 1;
    for n = 1:n_samples
        error(n) = y(n) - w(:,n)' * x(:,n);
        if n < length(x)
            w(:, n+1) = (1-mu*leak)*w(:, n) + mu*conj(error(n))*x(:,n);
        end
    end
end