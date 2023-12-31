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

% load ECG data & MR trigger markers
ecg = data(:, 1);
mrk = data(:, 4);

% Sampling rate in Hz
fs = 5000;

% TR in seconds
TR = 1.5;

% %% Parse ECG by MR marker
% 
% % Find the indices of positive and negative values in 'mrk'
% positive_indices = find(mrk > 0);
% negative_indices = find(mrk <= 0);
% 
% % Initialize a cell array to store the series of positive values
% positive_series = cell(0);
% 
% % Find the indices where positive series start
% start_indices = positive_indices([true; diff(positive_indices) > 1]);
% 
% % Iterate through the start indices to extract positive series
% for i = 1:length(start_indices)
%     start_index = start_indices(i);
% 
%     % Find the next negative index or the end of positive indices
%     if i < length(start_indices)
%         next_negative_index = find(negative_indices > start_indices(i), 1);
%         end_index = negative_indices(next_negative_index) - 1;
%     else
%         end_index = length(mrk);
%     end
% 
%     % Check if the length of the negative series is greater than 2 TR
%     if end_index - start_index + 1 > round(2*TR*fs)
%         positive_series{end+1} = start_index:end_index;
%     end
% end
% 
% % Extract the corresponding portions of the 'ecg' signal for each positive series
% ecg_localizer_series = cell(0);
% for i = 1:length(positive_series)
%     series_indices = positive_series{i};
%     ecg_localizer_series{i} = ecg(series_indices);
% end
% 
% % Combine the series into a single array
% ecg_localizer = cat(1, ecg_localizer_series{:});
% 
% % Create a time vector in seconds
% time = (0:length(ecg_localizer)-1) / fs;
% 
% % Plot the ecg_localizer_series
% plot(time, ecg_localizer)
% xlabel('Time (s)')
% ylabel('ECG')
% title('ECG Localizer Series')
% 
% keyboard
% 
% ecg = ecg_localizer;

%% Visualize raw ECG signal and PSD

% visualize raw input signal
figure; 
% ax1 = subplot(1,4,1);
plot(ecg); 
title('Unfiltered ECG Signal Lead I');
xlabel('time in ms'); ylabel('voltage in mV');
hold all; 

% Compute the power spectral density of the ECG signal
[pxx, f] = pwelch(ecg, [], [], [], fs, 'onesided');
pmax = pwelch(ecg,[],[],[],fs,'maxhold');
pmin = pwelch(ecg,[],[],[],fs,'minhold');

% Find the maximum and minimum values of the PSD % think this is redundant
% with lines 47-48 above

[max_pxx, max_idx] = max(pxx);
[min_pxx, min_idx] = min(pxx);

% Plot the power spectral density
figure;
%plot(f, 10*log10(pxx)); % unsure of difference between 10*log10(pxx) and
%pow2db(pxx), believe that they are equivalent approaches
hold on
plot(f,pow2db(pxx))
hold on
plot(f,pow2db([pmax pmin]),':')
hold off
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
legend('pwelch','maxhold','minhold')
grid on;
title('Power Spectral Density of ECG Signal');

%% Define and apply Zero-phase Butterworth IIR Bandpass Filter 
% https://www.mathworks.com/matlabcentral/answers/493660-how-do-you-design-your-ecg-bandpass

% % High-pass filter cutoff frequency (Hz) (test from 0.5hz)
% high_cutoff = 0.04;
% 
% % Low-pass filter cutoff frequency (Hz) 
% low_cutoff = 30;
% 
% % Powerline frequency (Hz)
% powerline_freq = 60;
% 
% % Design the Butterworth filters
% [b_high, a_high] = butter(2, high_cutoff/(fs/2), 'high');
% [b_low, a_low] = butter(2, low_cutoff/(fs/2), 'low');
% 
% % Design the notch filter
% % Note that 'notch' requires Signal Processing Toolbox R2016b or later
% wo = powerline_freq/(fs/2);  
% bw = wo/35;
% [b_notch, a_notch] = iirnotch(wo, bw);
% 
% % Apply the high-pass filter
% ecg_high = filtfilt(b_high, a_high, ecg);
% 
% % Apply the notch filter
% ecg_notch = filtfilt(b_notch, a_notch, ecg_high);
% 
% % Apply the low-pass filter
% ecg_filtered = filtfilt(b_low, a_low, ecg_notch);

% %% Test filter as per http://www.ijstr.org/final-print/dec2019/Digital-Iir-Filters-For-Heart-Rate-Variability-A-Comparison-Between-Butterworth-And-Elliptic-Filters.pdf
% 
% % Define the filter specifications
% passband_freq = 0.8;   % Passband frequency in Hz
% stopband_freq = 20;    % Stopband frequency in Hz
% passband_ripple = 0;   % Passband ripple in dB
% stopband_attenuation = 80;   % Stopband attenuation in dB
% 
% % Normalize the frequencies
% normalized_passband_freq = passband_freq / (fs/2);
% normalized_stopband_freq = stopband_freq / (fs/2);
% 
% % Design the Butterworth bandpass filter
% [b, a] = butter(10, [normalized_passband_freq, normalized_stopband_freq], 'bandpass');
% 
% % Apply the filter to the ECG signal
% ecg_filtered = filtfilt(b, a, ecg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test new butterworth to avoid producing negative voltages

% Define the filter parameters
high_cutoff = 0.04;  % High-pass filter cutoff frequency in Hz
low_cutoff = 50;   % Low-pass filter cutoff frequency in Hz

% Design the Butterworth filters
[b_high, a_high] = butter(2, high_cutoff/(fs/2), 'high');
[b_low, a_low] = butter(2, low_cutoff/(fs/2), 'low');

% Apply the high-pass filter
ecg_high = filtfilt(b_high, a_high, ecg);

% Apply the low-pass filter
ecg_filtered = filtfilt(b_low, a_low, ecg_high);

% Adjust the baseline voltage
ecg_baseline = min(ecg_filtered);
ecg_filtered_adjusted = ecg_filtered - ecg_baseline;
ecg_filtered = ecg_filtered_adjusted;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Test new filter - Kirstin filter
% 
% d = designfilt("lowpassfir", ...
%     PassbandFrequency=0.15,StopbandFrequency=0.2, ...
%     PassbandRipple=1,StopbandAttenuation=60, ...
%     DesignMethod="equiripple");
% ecg_filtered_low = filtfilt(d,ecg);
% %y1 = filter(ecg);
% 
% % Butterworth
% d1 = designfilt("lowpassiir",FilterOrder=12, ...
%     HalfPowerFrequency=0.15,DesignMethod="butter");
% ecg_filtered = filtfilt(d1,ecg_filtered_low);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the adjusted filtered ECG signal
figure;
plot(ecg_filtered);
title('Filtered ECG Signal with Adjusted Baseline');
xlabel('Time (s)');
ylabel('Voltage (mV)');

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
pmax_filt = pwelch(ecg_filtered,[],[],[],fs,'maxhold');
pmin_filt = pwelch(ecg_filtered,[],[],[],fs,'minhold');

% Plot the power spectral density of the filtered ECG signal
figure;
%plot(f_filt, 10*log10(pxx_filt));
hold on
plot(f_filt,pow2db(pxx_filt))
hold on
plot(f_filt,pow2db([pmax_filt pmin_filt]),':')
hold off
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
% [pks, locs] = findpeaks(ecg_filtered, 'MinPeakHeight', min(ecg_filtered) + 0.5 * (max(ecg_filtered) - min(ecg_filtered)), 'MinPeakDistance', 0.5 * fs);

[pks, locs] = findpeaks(ecg_filtered, 'MinPeakHeight', min(ecg_filtered) + 0.5 * (max(ecg_filtered) - min(ecg_filtered)), 'MinPeakDistance', 0.5 * fs, 'Annotate', 'extents');

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

%% Compute R-R interval tachogram (values not correct)

% Compute the time differences between consecutive R peaks (R-R intervals) in milliseconds
% rr_intervals = diff(locs) * 1000 / fs;
%rr_intervals = diff(x_locs);

peak_locs_ms = (locs - 1) * 1000 / fs;
rr_intervals = diff(peak_locs_ms);

% Generate time vector for x-axis
time = (locs(1:end-1) - 1) / fs;

% Plot the R-R interval tachogram
figure;
plot(time, rr_intervals, 'b.-');
title('R-R Interval Tachogram');
xlabel('Time (s)');
ylabel('R-R Interval Duration (ms)');

%% Filter RR Interval time series

% Apply filtering to remove implausible values
rr_intervals_filtered = rr_intervals;
%rr_intervals_filtered(rr_intervals < 500 | rr_intervals > 1800) = NaN; %
%double check this for 1500 vs 1800
rr_intervals_filtered(rr_intervals < 500 | rr_intervals > 1500) = NaN;

% Interpolate over removed values using cubic spline interpolation
rr_intervals_corrected = interp1(time, rr_intervals_filtered, time, 'pchip');

% %% Test tachogram band pass filter
% 
% % Apply the high-pass filter
% rr_intervals_corrected_high = filtfilt(b_high, a_high, rr_intervals_corrected);
% 
% % Apply the low-pass filter
% rr_intervals_corrected = filtfilt(b_low, a_low, rr_intervals_corrected_high);

% Plot the corrected R-R interval tachogram
figure;
plot(time, rr_intervals_corrected, 'b.-');
title('Corrected R-R Interval Tachogram');
xlabel('Time (s)');
ylabel('R-R Interval Duration (ms)');

%% Rerun filtering omitting noise sections

% Find the indices of non-NaN values in the RR intervals
valid_indices = find(~isnan(rr_intervals_corrected));

% Exclude RR intervals greater than 1500 ms
rr_intervals_corrected(rr_intervals_corrected > 1500) = NaN;

% Extract the valid time and RR intervals
time_valid = time(valid_indices);
rr_intervals_valid = rr_intervals_corrected(valid_indices);

% Remove NaN values from RR intervals
rr_intervals_valid = rr_intervals_valid(~isnan(rr_intervals_valid));

% Plot the corrected and valid R-R interval tachogram
figure;
plot(rr_intervals_valid, 'b.-');
title('Corrected and Valid R-R Interval Tachogram');
xlabel('Time (s)');
ylabel('R-R Interval Duration (ms)');

%% Time-domain HRV analysis

% Calculate average heart rate (beats per minute)
heart_rate = 60 / (mean(rr_intervals_valid) / 1000);

% Calculate standard deviation of RR intervals (SDNN)
sdnn = std(rr_intervals_valid);

% Calculate root mean square of successive RR interval differences (RMSSD)
rmssd = sqrt(mean(diff(rr_intervals_valid).^2));

% Calculate percentage of successive RR intervals that differ by more than 50 ms (pNN50)
pnn50 = sum(abs(diff(rr_intervals_valid)) > 50) / length(rr_intervals_valid) * 100;

% Display HRV results
disp('Heart Rate Variability Analysis Results:');
disp(['Average Heart Rate: ', num2str(heart_rate), ' beats per minute']);
disp(['SDNN: ', num2str(sdnn), ' milliseconds']);
disp(['RMSSD: ', num2str(rmssd), ' milliseconds']);
disp(['pNN50: ', num2str(pnn50), ' %']);

%% Frequency-domain HRV analysis

    % Define the desired frequency resolution
    frequency_resolution = 0.001; % Specify the desired frequency resolution in Hz

    % Define the frequency range of interest
    frequency_range = [.04, .4]; % Specify the frequency range of interest in Hz

    % Generate a higher resolution frequency vector
    freq_points = frequency_range(1):frequency_resolution:frequency_range(2);

% Perform power spectral density calculation only if valid RR intervals exist
if ~isempty(rr_intervals_valid)
    
    % Compute power spectral density (PSD) of RR intervals with increased resolution
    [psd, freq] = pwelch(rr_intervals_valid, [], [], freq_points, fs);
    pmax_psd = pwelch(rr_intervals_valid,[],[],freq_points,fs,'maxhold');
    pmin_psd = pwelch(rr_intervals_valid,[],[],freq_points,fs,'minhold');
    
    % Define frequency bands of interest
    lf_band = freq >= 0.04 & freq <= 0.15;
    hf_band = freq > 0.15 & freq <= 0.4;

    % Calculate LF power
    lf_power = trapz(psd(lf_band));

    % Calculate HF power
    hf_power = trapz(psd(hf_band));

    % Calculate total power
    total_power = trapz(psd);

    % Calculate LF/HF ratio
    lf_hf_ratio = lf_power / hf_power;

    % Calculate peak frequency for LF band
    [~, lf_peak_idx] = max(psd(lf_band));
    lf_peak_freq = freq(lf_band);
    lf_peak_freq = lf_peak_freq(lf_peak_idx);

    % Calculate peak frequency for HF band
    [~, hf_peak_idx] = max(psd(hf_band));
    hf_peak_freq = freq(hf_band);
    hf_peak_freq = hf_peak_freq(hf_peak_idx);

    % Calculate power in ms^2 for LF band
    lf_power_ms2 = lf_power * (1 / fs) * 1000;

    % Calculate power in ms^2 for HF band
    hf_power_ms2 = hf_power * (1 / fs) * 1000;

    % Calculate log power for LF and HF bands
    lf_log_power = 10 * log10(lf_power);
    hf_log_power = 10 * log10(hf_power);

    % desired_frequency = 0.04;  % Desired frequency in Hz
    % [~, index] = min(abs(freq - desired_frequency));  % Find the index of the closest frequency value

    % Calculate % total power for LF and HF bands
    lf_percent_total_power = (lf_power / total_power) * 100;
    hf_percent_total_power = (hf_power / total_power) * 100;

    % Calculate power in normalized units for LF and HF bands
    lf_normalized_power = lf_power / total_power;
    hf_normalized_power = hf_power / total_power;

    % Display PSD frequency-domain HRV results
    
    disp('PSD frequency domain analysis results.');
    disp(['LF Power: ', num2str(lf_power)]);
    disp(['LF Peak Frequency: ', num2str(lf_peak_freq), ' Hz']);
    disp(['LF Power (ms^2): ', num2str(lf_power_ms2), ' ms^2']);
    disp(['LF Log Power: ', num2str(lf_log_power), ' dB']);
    disp(['LF Power (% Total Power): ', num2str(lf_percent_total_power), ' %']);
    disp(['LF Power (Normalized Units): ', num2str(lf_normalized_power)]);
    
    disp(['HF Power: ', num2str(hf_power)]);
    disp(['HF Peak Frequency: ', num2str(hf_peak_freq), ' Hz']);
    disp(['HF Power (ms^2): ', num2str(hf_power_ms2), ' ms^2']);
    disp(['HF Log Power: ', num2str(hf_log_power), ' dB']);
    disp(['HF Power (% Total Power): ', num2str(hf_percent_total_power), ' %']);
    disp(['HF Power (Normalized Units): ', num2str(hf_normalized_power)]);

    % Plot the power spectral density (PSD) of RR intervals
    figure;
    %plot(freq, 10*log10(psd));
    hold on
    plot(freq,pow2db(psd))
    hold on
    %plot(freq,pow2db([pmax_psd pmin_psd])); %,':')
    hold off
    title('Power Spectral Density (PSD) of RR Intervals');
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB/Hz)');
    grid on;
else
    disp('No valid RR intervals for power spectral density calculation.');
end

%% TEST

% Perform FFT-based frequency domain analysis only if valid RR intervals exist
if ~isempty(rr_intervals_valid)
    
    % Compute the FFT of RR intervals
    fft_data = fft(rr_intervals_valid);
    
    % Compute the single-sided power spectrum
    power_spectrum = abs(fft_data).^2 / length(rr_intervals_valid);
    
    % Define the frequency range of interest
    min_freq = 0.04;  % Minimum frequency in Hz
    max_freq = 0.4;   % Maximum frequency in Hz
    
    % % Compute the corresponding frequency values
    % freq_resolution = fs / length(fft_data);  % Frequency resolution
    % freq = (0:length(fft_data)-1) * freq_resolution;

    % Define the desired frequency resolution
    %frequency_resolution = 0.001; % Specify the desired frequency resolution in Hz
    frequency_resolution = freq_points(2) - freq_points(1);

    % Generate a higher resolution frequency vector
    freq = (min_freq:frequency_resolution:max_freq);

    % Find the indices corresponding to the frequency range of interest
    freq_indices = freq >= min_freq & freq <= max_freq;
    
    % Find the indices corresponding to the frequency range of interest
    %freq_indices = freq >= min_freq & freq <= max_freq;
    
    % Extract the power spectrum and frequency values within the specified range
    power_spectrum_range = power_spectrum(freq_indices);
    freq_range = freq(freq_indices);
    
    % Calculate LF power
    lf_band = freq_range >= 0.04 & freq_range <= 0.15;
    lf_power = trapz(power_spectrum_range(lf_band));
    
    % Calculate HF power
    hf_band = freq_range > 0.15 & freq_range <= 0.4;
    hf_power = trapz(power_spectrum_range(hf_band));
    
    % Calculate total power
    total_power = trapz(power_spectrum_range);
    
    % Calculate LF/HF ratio
    lf_hf_ratio = lf_power / hf_power;
    
    % Calculate peak frequency for LF band
    [~, lf_peak_idx] = max(power_spectrum_range(lf_band));
    lf_peak_freq = freq_range(lf_band);
    lf_peak_freq = lf_peak_freq(lf_peak_idx);
    
    % Calculate peak frequency for HF band
    if any(hf_band)
        [~, hf_peak_idx] = max(power_spectrum_range(hf_band));
        hf_peak_freq = freq_range(hf_band);
        hf_peak_freq = hf_peak_freq(hf_peak_idx);
    else
        hf_peak_freq = NaN;
    end
    
    % Calculate power in ms^2 for LF band
    lf_power_ms2 = lf_power * (1 / fs) * 1000;
    
    % Calculate power in ms^2 for HF band
    hf_power_ms2 = hf_power * (1 / fs) * 1000;
    
    % Calculate log power for LF and HF bands
    lf_log_power = 10 * log10(lf_power);
    hf_log_power = 10 * log10(hf_power);
    
    % Calculate % total power for LF and HF bands
    lf_percent_total_power = (lf_power / total_power) * 100;
    hf_percent_total_power = (hf_power / total_power) * 100;
    
    % Calculate power in normalized units for LF and HF bands
    lf_normalized_power = lf_power / total_power;
    hf_normalized_power = hf_power / total_power;
    
    % Display FFT frequency-domain HRV results
    
    disp('FFT frequency domain analysis results.');
    disp(['LF Power: ', num2str(lf_power)]);
    disp(['LF Peak Frequency: ', num2str(lf_peak_freq), ' Hz']);
    disp(['LF Power (ms^2): ', num2str(lf_power_ms2), ' ms^2']);
    disp(['LF Log Power: ', num2str(lf_log_power), ' dB']);
    disp(['LF Power (% Total Power): ', num2str(lf_percent_total_power), ' %']);
    disp(['LF Power (Normalized Units): ', num2str(lf_normalized_power)]);
    
    disp(['HF Power: ', num2str(hf_power)]);
    disp(['HF Peak Frequency: ', num2str(hf_peak_freq), ' Hz']);
    disp(['HF Power (ms^2): ', num2str(hf_power_ms2), ' ms^2']);
    disp(['HF Log Power: ', num2str(hf_log_power), ' dB']);
    disp(['HF Power (% Total Power): ', num2str(hf_percent_total_power), ' %']);
    disp(['HF Power (Normalized Units): ', num2str(hf_normalized_power)]);
    
    % Plot the power spectrum
    figure;
    plot(freq_range, 10*log10(power_spectrum_range));
    title('Power Spectrum of RR Intervals');
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    grid on;
    
else
    disp('No valid RR intervals for frequency domain analysis.');
end

%% Test AR spectrum code

% Perform autoregressive (AR) modeling for frequency domain HRV analysis
if ~isempty(rr_intervals_valid)
    
    % Estimate the autoregressive (AR) model parameters
    order = 16; % Specify the order of the AR model
    ar_model = ar(rr_intervals_valid, order);
    
    % Compute the frequency response of the AR model
    freq_response = freqz(ar_model.a, ar_model.c, freq_points, fs);
    
    % Compute the power spectral density (PSD) from the AR model frequency response
    psd = abs(freq_response).^2;
    freq = freq_points;
    pmax_psd = max(psd);
    pmin_psd = min(psd);
    
    % Define frequency bands of interest
    lf_band = freq >= 0.04 & freq <= 0.15;
    hf_band = freq > 0.15 & freq <= 0.4;
    
    % Calculate LF power
    lf_power = trapz(psd(lf_band));

    % Calculate HF power
    hf_power = trapz(psd(hf_band));

    % Calculate total power
    total_power = trapz(psd);

    % Calculate LF/HF ratio
    lf_hf_ratio = lf_power / hf_power;

    % Calculate peak frequency for LF band
    [~, lf_peak_idx] = max(psd(lf_band));
    lf_peak_freq = freq(lf_band);
    lf_peak_freq = lf_peak_freq(lf_peak_idx);

    % Calculate peak frequency for HF band
    [~, hf_peak_idx] = max(psd(hf_band));
    hf_peak_freq = freq(hf_band);
    hf_peak_freq = hf_peak_freq(hf_peak_idx);

    % Calculate power in ms^2 for LF band
    lf_power_ms2 = lf_power * (1 / fs) * 1000;

    % Calculate power in ms^2 for HF band
    hf_power_ms2 = hf_power * (1 / fs) * 1000;

    % Calculate log power for LF and HF bands
    lf_log_power = 10 * log10(lf_power);
    hf_log_power = 10 * log10(hf_power);

    % Calculate % total power for LF and HF bands
    lf_percent_total_power = (lf_power / total_power) * 100;
    hf_percent_total_power = (hf_power / total_power) * 100;

    % Calculate power in normalized units for LF and HF bands
    lf_normalized_power = lf_power / total_power;
    hf_normalized_power = hf_power / total_power;

    % Display AR modeling frequency-domain HRV results
    
    disp('AR modeling frequency domain analysis results.');
    disp(['LF Power: ', num2str(lf_power)]);
    disp(['LF Peak Frequency: ', num2str(lf_peak_freq), ' Hz']);
    disp(['LF Power (ms^2): ', num2str(lf_power_ms2), ' ms^2']);
    disp(['LF Log Power: ', num2str(lf_log_power), ' dB']);
    disp(['LF Power (% Total Power): ', num2str(lf_percent_total_power), ' %']);
    disp(['LF Power (Normalized Units): ', num2str(lf_normalized_power)]);
    
    disp(['HF Power: ', num2str(hf_power)]);
    disp(['HF Peak Frequency: ', num2str(hf_peak_freq), ' Hz']);
    disp(['HF Power (ms^2): ', num2str(hf_power_ms2), ' ms^2']);
    disp(['HF Log Power: ', num2str(hf_log_power), ' dB']);
    disp(['HF Power (% Total Power): ', num2str(hf_percent_total_power), ' %']);
    disp(['HF Power (Normalized Units): ', num2str(hf_normalized_power)]);

    % Plot the power spectral density (PSD) of RR intervals
    figure;
    plot(freq, pow2db(psd));
    hold on
    %plot(freq, pow2db([pmax_psd, pmin_psd]), ':');
    hold off
    title('Power Spectral Density (PSD) of RR Intervals (AR Modeling)');
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB/Hz)');
    grid on;
else
    disp('No valid RR intervals for AR modeling and power spectral density calculation.');
end


%% Signal for saving

% Save signals into a new .mat file (update to draw filename from loaded ACQ .mat export file)
%save(fullfile(pathname, 'pilot003_run1_ECG.mat'), 'ecg', 'ecg_filtered', 'rr_intervals_valid');
