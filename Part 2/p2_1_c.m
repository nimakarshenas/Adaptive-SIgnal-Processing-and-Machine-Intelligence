clc; clear all; close all;
% realisation length
N = 100000;
% no. realisations
M = 100;
var = 0.25;
% get realisation
data = filter([1], [1, -0.1, -0.8], sqrt(var)*randn(N+500, M));
data = data(501:end, :); 
errors = zeros(2, M, N); 
step_sizes = [0.05, 0.01];
for i = 1:M
    for j = 1:2
        [~, errors(j, i, :)] = lms_estimator(data(:, i), 2, step_sizes(j)); 
    end
end
errors = squeeze(mean(errors(:, 100, 1000:end).^2, [2, 3]));
M1 = (errors(1) - var)/var;
M2 = (errors(2) - var)/var;


%% function for LMS estimation
function [params, error] = lms_estimator(data, order, step_size)

    params = zeros(order, length(data));
    error = zeros(size(data)) ;

    for i = order+1:length(data)
        current_error = data(i) - dot(flip(data(i-order:i-1)), params(:, i-1));
        error(i) = current_error;
        params(:, i) = params(:, i-1) + step_size*(current_error)*flip(data(i-order:i-1));
    end
end
