clc; clear all; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


complex_wind_data = zeros(3,5000,'like',1i);
files = {'high-wind','medium-wind','low-wind'};
cols=["#7E2F8E","#77AC30","#4DBEEE"];
lims = [4 2 0.5];
rhos = zeros(1, 3);
%% Plot wind data on complex plane
figure(1);
for i = 1:3
    data = load(files{i});
    complex_wind_data(i,:) = complex(data.v_east,data.v_north);
    str = split(files{i},'-');
    tit = strcat(upper(str{1}(1)), str{1}(2:end), ' Wind');
    subplot(1,3,i); 
    hold on; 
    set(gca,'fontsize', 16);
    scatter(real(complex_wind_data(i,:)), imag(complex_wind_data(i,:)), 3, ...
        'filled', 'MarkerFaceColor',cols(i), 'MarkerEdgeColor',cols(i))
    xlabel('$\mathcal{R}$'); ylabel('$\mathcal{I}$');
    grid on; grid minor;
    xlim([-lims(i) lims(i)]); ylim([-lims(i) lims(i)])
    title(tit);
    hold off;
    [rhos(i), ~] = circularity(complex_wind_data(i, :));
end

n_orders = 20;
MSPE = zeros(2,3,n_orders); 
step_sizes = [0.001, 0.007,0.1];
counter = 0;
figure(2);
for j = 1:3
    counter = counter + 1 ;
    subplot(1, 3,counter); hold on; set(gca,'fontsize', 18);
    for i = 1:2
        input = delayseq(complex_wind_data(j,:).',1);
        for l = 1:n_orders
            if i == 1
                [~,err,~] = clms(complex_wind_data(j,:).',input,l-1, step_sizes(j), 0);
            else
                [~,err,~] = aclms(complex_wind_data(j,:).',input,l-1, step_sizes(j), 0);
            end
            sq_err = abs(err).^2;
            MSPE(i,j,l) = mean(sq_err);
        end
        plot(1:n_orders, 10*log10(squeeze(MSPE(i,j,:))), LineWidth=1.5);
    end
    str = split(files{j},'-');
    set(gca, 'FontSize', 14)
    tit = strcat(upper(str{1}(1)), str{1}(2:end), ' Wind, $\mu=$', num2str(step_sizes(j)));
    legend('CLMS','ACLMS');
    xlabel('Filter Order'); ylabel('MSPE (dB)');
    grid on; grid minor;
    title(tit);
    hold off;
end


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