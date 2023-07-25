% Load EDA data
load('eda_test1.mat'); % Replace with your own data file

% Define filter parameters
fs = 5000; % Sampling frequency
f0 = 0.6667; % Notch frequency (in Hz)
f0_norm = f0/fs; % Notch frequency (normalized)
Q = 0.707; % Quality factor

% Calculate filter coefficients
w0 = 2*pi*f0_norm;
alpha = sin(w0)/(2*Q);
b0 = 1;
b1 = -2*cos(w0);
b2 = 1;
a0 = 1 + alpha;
a1 = -2*cos(w0);
a2 = 1 - alpha;

% Normalize coefficients by a0
b0 = b0/a0;
b1 = b1/a0;
b2 = b2/a0;
a1 = a1/a0;
a2 = a2/a0;
a0 = a0/a0;

% Create filter object
notchFilter = dsp.IIRFilter('Numerator', [b0 b1 b2], 'Denominator', [a0 a1 a2]);

% Apply filter to EDA data
eda_filtered = notchFilter(eda);

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