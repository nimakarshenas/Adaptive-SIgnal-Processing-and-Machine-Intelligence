%% Part a
close all; clear all; clc;


ecg1_og = readmatrix("Amir\s1.csv");
ecg1_og = ecg1_og(:, 3);
ecg2_og = readmatrix("Amir\s2.csv");
ecg2_og = ecg2_og(:, 3);

% ecg sampling freq
f_s = 500;


% plot to find artifact
figure(1);
subplot(1, 2, 1)
plot(ecg1_og(:, 1));
subplot(1, 2, 2)
plot(ecg2_og(:, 1));

% trim signals for each experiment
ecg1 = ecg1_og(4e3:222e3, 1);
ecg2 = ecg1_og(266e3:400e3, 1);
ecg3 = ecg2_og(5e3:142e3, 1);
ecg4 = ecg2_og(161e3:244e3, 1);
ecg5 = ecg2_og(266e3:366e3, 1);

[xRRI1, fs_RRI] = ECG_to_RRI_script(ecg1, f_s);
xRRI2 = ECG_to_RRI_script(ecg2, f_s);
xRRI3 = ECG_to_RRI_script(ecg3, f_s);
xRRI4 = ECG_to_RRI_script(ecg4, f_s);
xRRI5 = ECG_to_RRI_script(ecg5, f_s);


save('RRI_data', "fs_RRI", "xRRI1", "xRRI2", "xRRI3", "xRRI4", "xRRI5")