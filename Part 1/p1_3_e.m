close all
clear all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

model_order = 14;
n = 0:20;
nfft = 256;
N_repeat = 1000;
music_spectrum = zeros(N_repeat,nfft);


figure(1);
subplot(1,2,1)
for i=1:N_repeat
    
    noise = 0.05/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;

    [X,R] = corrmtx(x,model_order,'mod');
    [S,F] = pmusic(R,2,nfft,1,'corr');
    
    music_spectrum(i,:)=S;
    plot(F,S,'c'); set(gca,'xlim',[0.28 0.34]);
    hold on
end

mean_spectrum = mean(music_spectrum);
plot(F,mean_spectrum,'Linewidth',1);
grid on;
grid minor
xlabel('Normalised frequency ($\pi$ rad/sample)','FontSize',11); 
ylabel('PSD','FontSize',11);
title('\textbf{Realisations and Mean of MUSIC PSD estimate}','FontSize',11)

std_val = std(music_spectrum);
figure(1);
subplot(1,2,2)
plot(F,std_val,'Linewidth',1);
grid on; 
grid minor
xlabel('Normalised frequency ($\pi$ rad/sample)','FontSize',11); 
ylabel('$\sigma_n$ (PSD)','FontSize',11);
set(gca,'xlim',[0.28 0.34]);
title('\textbf{Standard Deviation of MUSIC PSD estimate}','FontSize',11)