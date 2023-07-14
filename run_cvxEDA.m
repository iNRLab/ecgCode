%function run_cvxEDA()

%% Clear workspace
clear all
close all

%% Load EDA data file
% Open the GUI for file selection - select eda file)
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

y = eda;
yn = zscore(y);
Fs = 5000;
[r, p, t, l, d, e, obj] = cvxEDA(yn, 1/Fs);
figure, hold all
tm = (1:length(y))'/Fs;
plot(tm, yn)
plot(tm, r)
plot(tm, p)
plot(tm, t)