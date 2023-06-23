%function loadECGfromACQ()

%% Clear workspace
clear all
close all

%% Load ECG data file
% Open the GUI for file selection - select ACQ exported .mat file with physio data)
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

% Find the maximum and minimum values of the PSD
[max_pxx, max_idx] = max(pxx);
[min_pxx, min_idx] = min(pxx);

% Plot the power spectral density
figure;
plot(f, 10*log10(pxx));
grid on;
title('Power Spectral Density of ECG Signal');
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');

%% Define and apply Zero-phase Butterworth IIR Bandpass Filter 
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

%% Compute the power spectral density of the filtered ECG signal
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

%% Label R peaks in filtered ECG (todo still not optimized)

% Find R peaks in the filtered ECG signal
[pks, locs] = findpeaks(ecg_filtered, 'MinPeakHeight', min(ecg_filtered) + 0.5 * (max(ecg_filtered) - min(ecg_filtered)), 'MinPeakDistance', 0.5 * fs);

% Compute the correct x-axis values for the R peaks
x_locs = (locs - 1) / fs;

% Plot the filtered ECG signal with labeled R peaks
figure;
plot((0:length(ecg_filtered)-1) / fs, ecg_filtered);
hold on;
plot(x_locs, pks, 'ro', 'MarkerSize', 8);
hold off;
title('Filtered ECG Signal with Labeled R Peaks');
xlabel('Time (s)');
ylabel('Voltage (mV)');
legend('Filtered ECG', 'R Peaks');

% Define the duration of the window in seconds
window_duration = 8;

% Compute the number of samples corresponding to the window duration
window_samples = window_duration * fs;

% Compute the end index of the window
end_idx = min(start_idx + window_samples - 1, length(ecg_filtered));

% Extract the filtered ECG signal within the specified window
window_ecg = ecg_filtered(start_idx:end_idx);

% Find the R peaks within the specified window
[pks_window, locs_window] = findpeaks(window_ecg, 'MinPeakHeight', min(window_ecg) + 0.5 * (max(window_ecg) - min(window_ecg)), 'MinPeakDistance', 0.5 * fs);

% Adjust the R peak locations to match the global time index
window_locs = locs_window + (start_idx - 1);

% Convert the adjusted R peak locations to the corresponding time values
window_time = t(window_locs) / 1000;  % Convert to seconds

% Plot the filtered ECG signal within the specified window
figure;
plot(t(start_idx:end_idx) / 1000, ecg_filtered(start_idx:end_idx));
hold on;
plot(window_time, pks_window, 'ro', 'MarkerSize', 8);
hold off;
title('Filtered ECG Signal with Labeled R Peaks');
xlabel('Time (s)');
ylabel('Voltage (mV)');
legend('Filtered ECG', 'R Peaks');

%% Heart Rate Variability Analysis (values not correct)

% Calculate R-R intervals using the global R peak locations
rr_intervals = diff(locs) / fs;

% Calculate average heart rate (beats per minute)
heart_rate = 60 / mean(rr_intervals);

% Perform additional HRV analysis (e.g., time-domain or frequency-domain analysis)

% Time-domain HRV analysis
% Calculate standard deviation of RR intervals (SDNN)
sdnn = std(rr_intervals);

% Calculate root mean square of successive RR interval differences (RMSSD)
rmssd = sqrt(mean(diff(rr_intervals).^2));

% Calculate percentage of successive RR intervals that differ by more than 50 ms (pNN50)
pnn50 = sum(abs(diff(rr_intervals)) > 0.05) / length(rr_intervals) * 100;

% Display HRV results
disp('Heart Rate Variability Analysis Results:');
disp(['Average Heart Rate: ', num2str(heart_rate), ' beats per minute']);
disp(['SDNN: ', num2str(sdnn), ' seconds']);
disp(['RMSSD: ', num2str(rmssd), ' seconds']);
disp(['pNN50: ', num2str(pnn50), ' %']);

% Frequency-domain HRV analysis
% Compute power spectral density (PSD) of RR intervals
[psd, freq] = pwelch(rr_intervals, [], [], [], fs);

% Extract frequency components of interest
lf_band = freq >= 0.04 & freq <= 0.15;
hf_band = freq > 0.15 & freq <= 0.4;

% Calculate LF power
lf_power = sum(psd(lf_band)) * (fs/length(psd));

% Calculate HF power
hf_power = sum(psd(hf_band)) * (fs/length(psd));

% Calculate total power
total_power = sum(psd) * (fs/length(psd));

% Calculate LF/HF ratio
lf_hf_ratio = lf_power / hf_power;

% Display frequency-domain HRV results
disp(['LF Power: ', num2str(lf_power)]);
disp(['HF Power: ', num2str(hf_power)]);
disp(['Total Power: ', num2str(total_power)]);
disp(['LF/HF Ratio: ', num2str(lf_hf_ratio)]);

%% Signal for saving

% Save signals into a new .mat file (update to draw filename from loaded ACQ .mat export file)
save(fullfile(pathname, 'pilot003_run1_ECG.mat'), 'ecg', 'ecg_filtered');
