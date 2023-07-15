clc; clear variables; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

n_samples = 1500;
f = [100*ones(500,1); ...
    100*ones(500,1) + ([501:1000]' - 500)/2; ...
    100*ones(500,1) + (([1001:1500]' -  1000)/25).^2];
phi = cumsum(f);
fs = 1500;
variance = 0.05;
eta = wgn(n_samples,1,pow2db(variance),'complex');
y = exp(1i*2*pi*phi/fs) + eta;
figure(1); subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot([1:n_samples], f, LineWidth=1.5);
ylim([0 500])
grid on; grid minor;
xlabel('Time (samples)'); ylabel('Frequency (Hz)');
title('Frequency of FM signal'); hold off;
order = [1, 5, 10, 20];


subplot(1,2,2); hold on; set(gca,'fontsize', 16);
H = zeros(4,1500); w = zeros(4,1500);
for i = 1:4
    a = aryule(y,order(i));
    [H(i,:),w(i,:)] = freqz(1,a,n_samples,fs);
    P = abs(H(i,:)).^2;
    plot(w(i,:), pow2db(P), LineWidth=1.5);
end
title('AR Spectrum estimates')
xlabel('Frequency (Hz)'); 
ylabel('PSD (dB)');
leg = arrayfun(@(a)strcat('p=',num2str(a)),order,'uni',0);
grid on; grid minor;
leg{end+1} = 'FMavg';
legend(leg); 
hold off;
