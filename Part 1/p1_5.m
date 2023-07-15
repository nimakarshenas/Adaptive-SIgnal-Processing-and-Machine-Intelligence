%% Part a
close all; clear all; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


load RRI_data.mat

xRRI1 = detrend(normalize(xRRI1'));
xRRI2 = detrend(normalize(xRRI2'));
xRRI3 = detrend(normalize(xRRI3'));

xRRI1 = detrend(normalize(xRRI1'));
xRRI2 = detrend(normalize(xRRI2'));
xRRI3 = detrend(normalize(xRRI3'));
xRRI4 = detrend(normalize(xRRI4'));
xRRI5 = detrend(normalize(xRRI4'));


RRI_data = {xRRI1; xRRI3; xRRI4};

figure(2);
for j = 1:length(RRI_data) % Trials
    subplot(1, length(RRI_data), j); hold on;
    L = length(RRI_data{j}); 
    wl = [L, 150*fs_RRI, 50*fs_RRI];
    colours = ["r", "b", "c"];
    for i = 1:3 % Window sizes
        [P, w] = pwelch(RRI_data{j}, wl(i), 0, L, fs_RRI);
        plot(w, 10*log10(P), LineWidth=i*0.5, Color=colours(i));
    end
    set(gca,'fontsize', 16);
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB)'); 
    grid on
    grid minor
    title(sprintf('Trial %d', j), FontSize=18);
    legend('Standard', '$\Delta_t$ = 150s', '$\Delta_t$ = 50s', FontSize=12); hold off;
end
sgtitle('\textbf{Standard and Averaged RRI PSD Estimates}', FontSize=20);

%% Part c

order = [6, 11, 12]; 

figure(3);
for i = 1:length(RRI_data)
    
    L = length(RRI_data{i});
    [P_og, w] = pwelch(RRI_data{i}, L, 0, L, fs_RRI);

    subplot(1, length(RRI_data), i)
    plot(w, 10*log10(P_og));
    hold on;
    [P_est, w] = pyulear(RRI_data{i}, order(i), 2048, fs_RRI); % power spectrum estimate given AR model
    plot(w, 10*log10(P_est), LineWidth=2);
 
    hold off;
    set(gca,'fontsize', 20);
    grid on;
    grid minor;
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB)'); 
    title(sprintf('Trial %d', i), FontSize=22);
    legend('Standard', sprintf('AR(%d)', order(i))); 
    hold off;
end
sgtitle("\textbf{RRI PSD Estimation via AR Modelling}", fontsize=24)
