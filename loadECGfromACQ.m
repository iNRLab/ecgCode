%function loadECGfromACQ()

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
disp(size(data))

% Plot ECG data
plot(data(:, 1)) 
ecg = data(:, 1);

% Sampling rate
fs = 5000;

% Signal for saving
signal = ecg;

% Save the signal into a new .mat file
save(fullfile(pathname, 'pilot003_ECG.mat'), 'signal');