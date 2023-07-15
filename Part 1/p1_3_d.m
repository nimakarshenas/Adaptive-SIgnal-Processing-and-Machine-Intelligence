close all
clear all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


f_s=1;
N=[20,40,80,100];

figure(1);
for i=1:4
    n = 0:N(i)-1;
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;
    K=N(i)/1024;
    [pxx,f] = periodogram(x,rectwin(N(i)),round(N(i)/K),f_s);
    plot(f,10*log(pxx),'Linewidth',1)
    hold on
end

xlabel('Frequency (Hz)','FontSize',11)
ylabel('PSD (dB)','FontSize',11)
grid on
grid minor 
title('\textbf{PSD estimates of complex exponentials}','FontSize',11)
legend('N=20','N=40','N=80','N=100','FontSize',9)
xlim([0.24 0.38])