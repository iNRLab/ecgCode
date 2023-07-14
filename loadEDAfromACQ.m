%function loadEDAfromACQ()

%% Clear workspace
clear all
close all

%% Load EDA data file
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

% load EDA data & MR trigger markers
eda = data(:, 3);
mrk = data(:, 4);

% Sampling rate in Hz
fs = 5000;

% TR in seconds (Note Resting State 2.0 seconds and Task 1.5 sec have different TRs)
TR = 2.0;

%% Visualize raw EDA signal and PSD

% visualize raw input signal
figure; 
% ax1 = subplot(1,4,1);
plot(eda); 
title('Unfiltered EDA Signal');
xlabel('time in ms'); ylabel('voltage in mS');
hold all; 

% Compute the power spectral density of the ECG signal
[pxx, f] = pwelch(eda, [], [], [], fs, 'onesided');
pmax = pwelch(eda,[],[],[],fs,'maxhold');
pmin = pwelch(eda,[],[],[],fs,'minhold');

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
title('Power Spectral Density of EDA Signal');

%% Signal for saving

%Save signals into a new .mat file (update to draw filename from loaded ACQ .mat export file)
save(fullfile(pathname, 'LRN001_rest_run1_EDA.mat'), 'eda');
