%function loadECGfromACQ()
clear all

% Open the GUI for file selection
[filename, pathname] = uigetfile('*.mat', 'Select a MAT-file', 'MultiSelect', 'on');

% Check if a file is selected
if isequal(filename, 0)
    disp('No file was selected.')
    return
end

% Check if the selected file exists
if ~exist(fullfile(pathname, filename), 'file')
    disp('Selected file does not exist.')
    return
end

% Load selected .mat file
load(fullfile(pathname, filename));
cd(pathname)

% Display the size of loaded data
disp(size(data));

% load ECG data
ecg = data(:, 1);

% visualize raw input signal
figure; 
% ax1 = subplot(1,4,1);
plot(ecg); 
title('Unfiltered ECG Signal Lead I');
xlabel('time in ms'); ylabel('voltage in mV');
hold all; 

% Sampling rate
fs = 5000;

% Compute the power spectral density of the ECG signal
[pxx, f] = pwelch(ecg, [], [], [], fs, 'onesided');

% Plot the power spectral density
figure;
plot(f, 10*log10(pxx));
grid on;
title('Power Spectral Density of ECG Signal');
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');

% High-pass filter cutoff frequency (Hz)
high_cutoff = 0.5;

% Low-pass filter cutoff frequency (Hz)
low_cutoff = 30;

% Powerline frequency (Hz)
powerline_freq = 60;

% Design the Butterworth filters
[b_high, a_high] = butter(2, high_cutoff/(fs/2), 'high');
[b_low, a_low] = butter(2, low_cutoff/(fs/2), 'low');

% Design the notch filter
% Note that 'notch' requires Signal Processing Toolbox R2016b or later
wo = powerline_freq/(fs/2);  
bw = wo/35;
[b_notch, a_notch] = iirnotch(wo, bw);

% Apply the high-pass filter
ecg_high = filtfilt(b_high, a_high, ecg);

% Apply the notch filter
ecg_notch = filtfilt(b_notch, a_notch, ecg_high);

% Apply the low-pass filter
ecg_filtered = filtfilt(b_low, a_low, ecg_notch);

% Plot filtered signal
figure; 
plot(ecg_filtered); 
title('Filtered ECG Signal Lead I');
xlabel('time in ms'); ylabel('voltage in mV');

% Plot original and filtered signal in the same figure
figure; 
hold on;
plot(ecg);
plot(ecg_filtered);
title('Original and Filtered ECG Signals');
xlabel('time in ms'); ylabel('voltage in mV');
legend('Unfiltered ECG', 'Filtered ECG');
hold off;

% Compute the power spectral density of the filtered ECG signal
[pxx_filt, f_filt] = pwelch(ecg_filtered, [], [], [], fs, 'onesided');

% Plot the power spectral density of the filtered ECG signal
figure;
plot(f_filt, 10*log10(pxx_filt));
grid on;
title('Power Spectral Density of Filtered ECG Signal');
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');

% Define time vector
t = (0:length(ecg)-1) / fs * 1000;  % Convert to milliseconds

% Choose a starting point for the plot (in ms)
start_time = 5000;  % For example

% Find the corresponding index for the start time
start_idx = find(t >= start_time, 1);

% Compute the index for the end of the 2000ms window
end_idx = start_idx + 2000 * fs / 1000 - 1;

% Make sure end_idx doesn't exceed the length of the signal
end_idx = min(end_idx, length(ecg));

% Plot original ECG signal over 2000ms window
figure;
plot(t(start_idx:end_idx), ecg(start_idx:end_idx));
title('Original ECG Signal over 2000ms Window');
xlabel('Time (ms)');
ylabel('Voltage (mV)');

% Plot filtered ECG signal over 2000ms window
figure;
plot(t(start_idx:end_idx), ecg_filtered(start_idx:end_idx));
title('Filtered ECG Signal over 2000ms Window');
xlabel('Time (ms)');
ylabel('Voltage (mV)');

% Create a new figure
figure;

% Plot original ECG signal over 2000ms window
plot(t(start_idx:end_idx), ecg(start_idx:end_idx));
hold on;

% Plot filtered ECG signal over 2000ms window
plot(t(start_idx:end_idx), ecg_filtered(start_idx:end_idx));

% Add title, labels, and legend
title('Original and Filtered ECG Signals over 2000ms Window');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
legend('Original', 'Filtered');

% Signal for saving

% Save signals into a new .mat file
save(fullfile(pathname, 'pilot003_run1_ECG.mat'), 'ecg', 'ecg_filtered');