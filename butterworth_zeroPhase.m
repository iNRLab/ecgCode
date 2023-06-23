% Sampling rate (Hz)
fs = 250;

% ECG data
% Assuming ecg is your ECG data
% ecg = your_data_here;

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