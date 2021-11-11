% The sampling rate is 1000 Hz
FS = 1000;

% Load ECG 1 into Nx1 vector from the file ecg_signal_1.dat
ecg1 = load("ecg_signal_1.dat");

% Load ECG 2 into Nx1 vector from the file ecg_signal_2.dat
ecg2 = load("ecg_signal_2.dat");

%% Power Spectrum

% Compute ECG 1 power spectrum
N1 = numel(ecg1);

fft_ecg1 = fft(ecg1,N1);

P_ecg1 = (fft_ecg1 .* conj(fft_ecg1))/N1;

% % Compute ECG 2 power spectrum
N2 = numel(ecg2);
fft_ecg2 = fft(ecg2,N2);
P_ecg2 = (fft_ecg2 .* conj(fft_ecg2)) / N2;

% % Compute power spectrum frequency bins from 0 Hz to the Nyquist frequency
f1 = 0:FS/numel(P_ecg1):FS/2;
f2 = 0:FS/numel(P_ecg2):FS/2;

%% Moving Average Filter

windowSize = 10;
b = (1/windowSize)*ones(1,windowSize);
a = 1;

% Do the filtering using a, b, and ecg1
% For ecg1
ecg1_filtered = filter(b,a,ecg1);
% ...and ecg2
ecg2_filtered = filter(b,a,ecg2);

%% Derivative Based Filter

% Create derivative based filter coefficients a and b:
b_ = FS*[1, -1];
a = [1, -0.995];

normalizing_term = max(freqz(b_,a,2000));

b = b_ * (1/normalizing_term);

% Do the filtering using a, b, and ecg1
% For ecg1
ecg1_filtered = filter(b,a,ecg1);
% ...and ecg2
ecg2_filtered = filter(b,a,ecg2);


%% Comb Filter

% Create comb filter coefficients a and b:
b = [0.6310, -0.2149, 0.1512, -0.1288, 0.1227, -0.1288, 0.1512, -0.2149, 0.6310];
a = 1;
ecg1_filtered = filter(b,a,ecg1);
ecg2_filtered = filter(b,a,ecg2);

%% Cascading Filter

windowSize = 10;
b1 = (1/windowSize)*ones(1,windowSize);
a1 = 1;

b2_ = FS*[1, -1];
a2 = [1, -0.995];

normalizing_term = max(freqz(b2_,a2,2000));

b2 = b2_ * (1/normalizing_term);

b3 = [0.6310, -0.2149, 0.1512, -0.1288, 0.1227, -0.1288, 0.1512, -0.2149, 0.6310];
a3 = 1;

a = conv(a1,conv(a2,a3));
b = conv(b1,conv(b2,b3));

ecg1_filtered = filter(b,a,ecg1);
ecg2_filtered = filter(b,a,ecg2);

%% Plot

plot(1:1:numel(ecg1), ecg1_filtered);
hold on
plot(1:1:numel(ecg1), ecg1);
