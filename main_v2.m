% Test the state-space multitaper spectral estimation algorithm on
% SED10.mat (EEG data on a human subject under general anesthesia)
%
%   This code implements the state-space multitaper spectrogram  
%   described in Kim et al., 2018 PNAS. 
%
%   Usage:
%   main.m: Main code    
%   EM_parameters.m: Compute noise & state variance using EM algorithm 
%   periodogram.m: Compute periodogram
%   multitaper.m: Compute multitaper spectrogram
%   SS_ST.m: Compute SS periodogram
%   SS_MT.m: Compute SS mutitaper spectrogram
%
%   From the paper:
%  "State-space multitpaer time-freqeuncy analysis"
%   Kim, S-E, Behr, MK, Ba, D & Brown, EN
%   PNAS, 2018
%
%   Copyright 2018 The General Hospital Coporation, authored by Seong-Eun Kim, Ph.D.
%   This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
%   (http://creativecommons.org/licenses/by-nc-sa/4.0/)
%
%   Original author: Sebastian Gallo
%   Last modified 04/08/2025 by Christian Guay
%   galloseb@mit.edu
%   christian.guay@hsc.utah.edu
%
%************************************************************************** 
%% 1. Setup Environment

% clear; close all; clc;

dataset_name = 'ACP_Concatenated.mat';  % File name of the .mat dataset to be loaded. Use the flattenConcatEDG.m helper function to concatenate the multiple EDF's from a Sedline Root download first, and save as .mat

% Define experiment start time (HH:MM:SS)
experiment_start = [07, 49, 08];  % 08:49:22 (HH, MM, SS)

% Define desired start and end times (HH:MM:SS)
desired_start = [07, 49, 50];   % 10:06:00
desired_end = [12, 18, 44];    % 11:36:00

% Setting variable parameters
channel = 1; % Electrode we're using
fs = 178;    % Sampling frequency (Hz). Sedline fs = 178
fmax = 20;   % Max freq to analyze
cmin = -10;  % Min value in dB for spectral analysis
cmax = 15;   % Max value in dB for spectral analysis
win_length = 2; % length of window (second)
Individual_Spectrogram = true; %make false if want all spectrograms/periodograms

% Setting AMI and SMI parameters
threshold_type = 'baseline'; %baseline (50th percentile of set baseline window) vs. absolute (mcv)
baseline_start = 300; %time (s) where the baseline starts for automated AMI/SMI threshold calculations
baseline_end = 600; %time (s) where the baseline ends for automated AMI/SMI threshold calculations
AMI_threshold = 5; %absolute threshold for AMI
SMI_threshold = 5; %absolute threshold for SMI

% Channel column in MATLAB - Sedline channel - 10-20 Channel
% 1 - R2 - F8
% 2 - R1 - Fp2
% 3 - L1 - Fp1
% 4 - L2 - F7

%% Converting times and loading data

% Convert all times to seconds from experiment start (absolute time)
start_time = (desired_start(1) - experiment_start(1)) * 3600 + ...
             (desired_start(2) - experiment_start(2)) * 60 + ...
             (desired_start(3) - experiment_start(3));

end_time = (desired_end(1) - experiment_start(1)) * 3600 + ...
           (desired_end(2) - experiment_start(2)) * 60 + ...
           (desired_end(3) - experiment_start(3));

disp('Initializing State-Space Multitaper Spectrogram Analysis...');

data_path = fullfile('data', dataset_name);

% Load EEG data 
data = load_eeg_data(data_path);

% Load EEG data
disp(['Loading EEG data from ', dataset_name, '...']);


% Extract EEG data
    % eeg_data = data(:, elec);
    start_idx = round(start_time * fs) + 1;  % Convert and adjust to 1-based index
    end_idx = round(end_time * fs) + 1;
    eeg_data = data(start_idx:end_idx, channel);
    num_samples = length(data);  % Total number of EEG samples

    % Generate time vector in seconds
    time_seconds = (0:num_samples-1) / fs;
    
    %Generate time vector in minutes
    time_minutes = time_seconds/60;
    
    % Trim the time vector to match EEG data length
    trimmed_time = time_minutes(start_idx:end_idx);


disp('EEG data loaded successfully.');


%% 2. Preprocessing
disp('Preprocessing EEG data...');
Nt = length(eeg_data);
sf = 1/win_length; % one step of frequency
nw = win_length*fs; % the number of elements in a window
N = floor(Nt/nw); % the number of window
disp(['number of windows:', num2str(N)]);

% Matrix form of data according to the size of window
yy = reshape(eeg_data(1:nw*N),nw,N); 
rtaper = rectwin(nw);
rtaper = rtaper/sqrt(nw);

kyy = zeros(nw,N);
for i = 1 : N
    kyy(:,i) = rtaper.*yy(:,i);
end

% Fourier transfrom of data
frequencyY = fft(kyy,nw,1);

disp('EEG data pre-processed successfully.');
%% 3. EM Algorithm, reduces bleeding from other frequencies
disp('Performing EM algorithm data...');
% First we can limit the frequency range to 0 to (desired) Hz and we can adjust
% the max level depending on the EEG data for greater denoising. 
OBSNOISE_CUTOFF = 30*win_length; % 30 Hz

% Initially we can set the alpha and beta as 1 
alpha = 1;
beta = 1;

% Initial guess for the observation noise
observationNoise = 100;
GUESS_WINDOW_LENGTH = 150; % EM estimation: 5 min = 300 sec = 300/win = 150
if GUESS_WINDOW_LENGTH > size(frequencyY, 2)
    GUESS_WINDOW_LENGTH = size(frequencyY, 2);
    disp('Guess window (150) greater than actual');
end

% Estimation of parameters using the EM algorithm for non-tapered data
[sn, on, is, iv, lls] = EM_parameters(alpha, beta, ...
    frequencyY(:,1:GUESS_WINDOW_LENGTH), ...
    observationNoise, 1e-5, OBSNOISE_CUTOFF, 1000);
                                  
% Multitapering
TW = 1; % Time-bandwidth production
K = 3; % The Number of tapers
[tapers,concentrations]=dpss(nw,TW,K); % Get the optimal tapers

y_ex = repmat(yy,1,1,K);
y_ex = permute(y_ex,[1 3 2]);
mtY = y_ex;
for i = 1:N
    mtY(:,:,i) = tapers.*y_ex(:,:,i);
end
mtFrequencyY = fft(mtY,nw,1);

% Estimation of parameters using the EM algorithm for multtitapered data
[mtSn, mtOn, mtIs, mtIv, mtLls] = EM_parameters(alpha, beta, ...
    mtFrequencyY(:,:,1:GUESS_WINDOW_LENGTH), ...
    observationNoise, 1e-5, OBSNOISE_CUTOFF, 1000);
%% 4. Spectral Estimates
% periodogram
spect1 = periodogram(yy, fs);
% multitaper spectrogram
[spect2,spect2_taper, spect2_fc] = multitaper_fc(yy, fs, TW, K);
% state-space periodogram
spect3 = SS_ST(yy, fs, sn, on, is, iv);
% state-space multitaper spectrogram
[spect4, results_MT] = SS_MT_cov_v3(yy, fs, TW, K, mtSn, mtOn, mtIs, mtIv);

disp('Spectral estimates computed.');

%% 5. Plot Spectrograms
disp('Plotting results...');

if Individual_Spectrogram
        spectrograms = {...
        {'Multitaper Spectrogram', spect2}
    };  
else 
        spectrograms = {...
        {'Periodogram', spect1}, ...
        {'Multitaper Spectrogram', spect2}, ...
        {'State-Space periodogram', spect3}, ...
        {'State-Space Multitaper Spectrogram', spect4}
    };  
end      

if exist('trimmed_time', 'var')
    time_vector = trimmed_time;
else
    time_vector = [];
end

plot_spectrograms(spectrograms, N, win_length, fs, Individual_Spectrogram,time_vector,fmax,cmin,cmax);

%% 6. Plot raw EEG

figure
title('Raw EEG');
axis xy;
ylabel('Amplitude (mcv)');
xlabel ('Time (min)');
plot (trimmed_time, eeg_data)
ylim([-50 50]); % Adjusted as needed


%% 7. Calculate and Plot AMI and SMI

addpath(fullfile(matlabroot, 'toolbox', 'matlab', 'graphics'), '-begin');


% Convert baseline times to indices in the trimmed data
baseline_start_idx = baseline_start - start_time;
baseline_end_idx = baseline_end - start_time;

% Ensure baseline times are within bounds
if baseline_start_idx < 0
    baseline_start_idx = 0;
    disp('Warning: Adjusted baseline_start to beginning of trimmed data');
end
if baseline_end_idx > length(eeg_data)/fs
    baseline_end_idx = length(eeg_data)/fs;
    disp('Warning: Adjusted baseline_end to end of trimmed data');
end

if threshold_type == 'baseline'

%Calculate SMI and SMI with thresholds based on 50th percentile of baseline
    [AMI, SMI] = computeAMI_SMI(eeg_data, fs, baseline_start_idx, baseline_end_idx);
else

% Calculate AMI and SMI with absolute thresholds
[AMI, SMI] = computeAMI_SMI(eeg_data, fs, baseline_start_idx, baseline_end_idx, ...
    'AlphaThresholdType', 'absolute', 'AlphaThresholdValue', AMI_threshold, ...
    'SMIThresholdType', 'absolute', 'SMIThresholdValue', SMI_threshold);
end

% Plot AMI and SMI - using separate figures to avoid subplot issues
figure('Name', 'Alpha Modulation Index (AMI)');
plot(trimmed_time, AMI);
title('Alpha Modulation Index (AMI)');
xlabel('Time (min)');
ylabel('AMI');
ylim([0 1]);
grid on;

figure('Name', 'Slow Modulation Index (SMI)');
plot(trimmed_time, SMI);
title('Slow Modulation Index (SMI)');
xlabel('Time (min)');
ylabel('SMI');
ylim([0 1]);
grid on;

%% Function: Plot Spectrograms
function plot_spectrograms(spectrograms, N, win, fs, Individual_Spectrogram, time_vector,fmax,cmin,cmax)

    num_plots = length(spectrograms);  % Dynamically determine number of spectrograms
    
    % If time_vector is missing or empty, use default (1:N) * win / 60
    if nargin < 6 || isempty(time_vector)
        disp('No time vector provided. Using default time scale (1:N) * win / 60.');
        time_vector = (1:N) * win / 60;  % Default time in minutes
    end

    % Figure settings based on Individual_Spectrogram flag
    if Individual_Spectrogram
        fig_width = 1200;
        fig_height = 400;
        figure('Position', [100 100 fig_width fig_height]);
    else
        figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 0.7 0.9]);
    end

    % Loop through spectrograms and plot them
    for i = 1:num_plots
        subplot(num_plots, 1, i);  % Dynamically create subplots based on spectrogram count
        title_text = spectrograms{i}{1};  % Extract subplot title
        spect_data = spectrograms{i}{2};  % Extract spectrogram data

        % Plot using time_vector (either user-provided or default)
        imagesc(time_vector, (0:fmax*win)/win, pow2db(spect_data(1:fmax*win,:)));
        
        axis xy;
        set(gca, 'clim', [cmin cmax]);
        ylabel('Frequency (Hz)');
        colorbar;
        colormap jet;
        title(title_text);
    end
    
    xlabel('Time (minutes)');  % Time is always in minutes
    disp('Plotting complete.');
end

