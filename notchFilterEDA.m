% Load EDA data
load('eda_test1.mat'); % Replace with your own data file

% Define filter parameters
fs = 5000; % Sampling frequency
f0 = 0.667; % Notch frequency (in Hz)
Q = 50; % Quality factor

% Apply notch filter to EDA data
eda_filtered = notchFilterMRI(eda, fs, f0, Q);

% Plot original and filtered signals
t = (0:length(eda)-1)/fs;
figure;
subplot(2,1,1);
plot(t, eda);
title('Original EDA Signal');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2);
plot(t, eda_filtered);
title('Filtered EDA Signal');
xlabel('Time (s)');
ylabel('Amplitude');