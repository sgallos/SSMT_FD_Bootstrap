%   This code implements the state-space multitaper spectrogram  
%   described in Kim et al., 2018 PNAS, and AMI and SMI described in Adam 
%   et al., 2023 PNAS; and generates and interactive figure GUI to visualize
%   the data.  
%
%   Usage:
%   main.m: Main code    
%   EM_parameters.m: Compute noise & state variance using EM algorithm 
%   periodogram.m: Compute periodogram
%   multitaper.m: Compute multitaper spectrogram
%   SS_ST.m: Compute SS periodogram
%   SS_MT.m: Compute SS mutitaper spectrogram
%
%   References: 
%  "State-space multitpaer time-freqeuncy analysis"
%   Kim, S-E, Behr, MK, Ba, D & Brown, EN
%   PNAS, 2018
%
%   Copyright 2018 The General Hospital Coporation, authored by Seong-Eun Kim, Ph.D.
%   This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
%   (http://creativecommons.org/licenses/by-nc-sa/4.0/)
%
%   Original author: Sebastian Gallo
%   Last modified 04/28/2025 by Christian Guay
%   galloseb@mit.edu
%   christian.guay@hsc.utah.edu
%
%************************************************************************** 
%% 1. Setup Environment

clear; clc;

dataset_name = '20250509_TibialFracture.mat';  % File name of the .mat dataset to be loaded. Use the flattenConcatEDG.m helper function to concatenate the multiple EDF's from a Sedline Root download first, and save as .mat

% Define experiment start time (HH:MM:SS)
experiment_start = [07, 00, 00];  % 08:49:22 (HH, MM, SS)

% Define desired start and end times (HH:MM:SS)
desired_start = [07, 00, 00];   % 10:06:00
desired_end = [08, 22, 25];    % 11:36:00

% Add option for full file import
use_full_file = true; % Set to true to import entire file, false to use time ranges

% Setting variable parameters
channel = 1; % Electrode we're using
fs = 178;    % Sampling frequency (Hz). Sedline fs = 178, Ripple fs = whatever we downsample to
artifact_threshold = 100; % Threshold in microvolts for absolute amplitude rejection
fmax = 40;   % Max freq to analyze
cmin = -10;  % Min value in dB for spectral analysis
cmax = 15;   % Max value in dB for spectral analysis
win_length = 2; % length of window (second)
Individual_Spectrogram = true; %make false if want all spectrograms/periodograms

% Setting AMI and SMI parameters
threshold_type = 'baseline'; %baseline (50th percentile of set baseline window) vs. absolute (mcv)
baseline_start = 60;  %time (s) where the baseline starts (1 minute)
baseline_end = 120;   %time (s) where the baseline ends (2 minutes)
AMI_threshold = 5;    %absolute threshold for AMI (used only if threshold_type is 'absolute')
SMI_threshold = 5;    %absolute threshold for SMI (used only if threshold_type is 'absolute')

% Channel column in MATLAB - Sedline channel - 10-20 Channel
% 1 - R2 - F8
% 2 - R1 - Fp2
% 3 - L1 - Fp1
% 4 - L2 - F7

%% 2. Converting times and loading data


disp('Initializing time scaling and data loading.');

data_path = fullfile('data', dataset_name);

% Load EEG data
disp(['Loading EEG data from ', dataset_name, '...']);

if exist(data_path, 'file')
    loaded_data = load(data_path);
    % Get the first field name in the loaded structure
    fields = fieldnames(loaded_data);
    if ~isempty(fields)
        data = loaded_data.(fields{1});
        disp('Data loaded successfully.');
    else
        error('Loaded .mat file does not contain any data fields');
    end
else
    error(['File not found: ', data_path]);
end

% Handle data loading based on user preference
if use_full_file
    % Use the entire file
    disp('Loading entire EEG file...');
    eeg_data = data(:, channel);
    num_samples = length(data);
    
    % Generate time vector in seconds and minutes
    time_seconds = (0:num_samples-1) / fs;
    trimmed_time = time_seconds / 60; % Convert to minutes
else
    % Convert times to seconds from experiment start (absolute time)
    start_time = (desired_start(1) - experiment_start(1)) * 3600 + ...
                 (desired_start(2) - experiment_start(2)) * 60 + ...
                 (desired_start(3) - experiment_start(3));
    
    end_time = (desired_end(1) - experiment_start(1)) * 3600 + ...
               (desired_end(2) - experiment_start(2)) * 60 + ...
               (desired_end(3) - experiment_start(3));
    
    % Extract specific time range from EEG data
    start_idx = round(start_time * fs) + 1;  % Convert and adjust to 1-based index
    end_idx = round(end_time * fs) + 1;
    eeg_data = data(start_idx:end_idx, channel);
    num_samples = length(data);
    
    % Generate time vector in seconds
    time_seconds = (0:num_samples-1) / fs;
    
    % Generate time vector in minutes
    time_minutes = time_seconds/60;
    
    % Trim the time vector to match EEG data length
    trimmed_time = time_minutes(start_idx:end_idx);
end

disp('EEG data loaded successfully.');

%% 3. Preprocessing
disp('Preprocessing EEG data...');

% Store a copy of the original data for reference if needed
eeg_data_original = eeg_data;

% DC offset removal (mean subtraction)
eeg_data = eeg_data - mean(eeg_data);
disp('DC offset removed.');

% Simple artifact rejection based on amplitude thresholds
artifacts = abs(eeg_data) > artifact_threshold;
artifact_percentage = sum(artifacts) / length(eeg_data) * 100;
disp(['Detected ' num2str(artifact_percentage, '%.2f') '% of samples as artifacts.']);

% Replace artifacts with interpolation
if sum(artifacts) > 0
    artifact_indices = find(artifacts);
    valid_indices = find(~artifacts);
    eeg_data(artifact_indices) = interp1(valid_indices, eeg_data(valid_indices), artifact_indices, 'pchip');
    disp('Artifacts interpolated.');
end

% Continue with existing preprocessing
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


%% 4. EM Algorithm, reduces bleeding from other frequencies
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
%% 5. Spectral Estimates
% periodogram
spect1 = periodogram(yy, fs);
% multitaper spectrogram
[spect2,spect2_taper, spect2_fc] = multitaper_fc(yy, fs, TW, K);
% state-space periodogram
spect3 = SS_ST(yy, fs, sn, on, is, iv);
% state-space multitaper spectrogram
[spect4, results_MT] = SS_MT_cov_v3(yy, fs, TW, K, mtSn, mtOn, mtIs, mtIv);

disp('Spectral estimates computed.');

%% 6. Compute Aperiodic Component (1/f)
disp('Computing aperiodic component (1/f)...');

% Extract aperiodic component
[spectral_exponent, aperiodic_spectra, aperiodic_freqs] = extractAperiodicComponent(eeg_data, fs, win_length, 50);

% Convert spectral exponent time to match spectrogram time
exponent_time = linspace(trimmed_time(1), trimmed_time(end), length(spectral_exponent));

% Store for visualization
aperiodic_data.spectral_exponent = spectral_exponent;
aperiodic_data.exponent_time = exponent_time;
aperiodic_data.aperiodic_spectra = aperiodic_spectra;
aperiodic_data.freqs = aperiodic_freqs;

disp('Aperiodic component computed successfully.');

%% 7. Compute AMI and SMI

addpath(fullfile(matlabroot, 'toolbox', 'matlab', 'graphics'), '-begin');

% Convert baseline times to indices based on the import mode
if use_full_file
    % For full file mode, baseline_start and baseline_end are already in seconds
    baseline_start_idx = baseline_start;
    baseline_end_idx = baseline_end;
else
    % For time-constrained mode, adjust baseline times relative to the segment start
    baseline_start_idx = baseline_start - start_time;
    baseline_end_idx = baseline_end - start_time;
end

% Ensure baseline times are within bounds
if baseline_start_idx < 0
    baseline_start_idx = 0;
    disp('Warning: Adjusted baseline_start to beginning of trimmed data');
end
if baseline_end_idx > length(eeg_data)/fs
    baseline_end_idx = length(eeg_data)/fs;
    disp('Warning: Adjusted baseline_end to end of trimmed data');
end

if strcmp(threshold_type, 'baseline')
    % Calculate SMI and SMI with thresholds based on 50th percentile of baseline
    [AMI, SMI] = computeAMI_SMI(eeg_data, fs, baseline_start_idx, baseline_end_idx);
else
    % Calculate AMI and SMI with absolute thresholds
    [AMI, SMI] = computeAMI_SMI(eeg_data, fs, baseline_start_idx, baseline_end_idx, ...
        'AlphaThresholdType', 'absolute', 'AlphaThresholdValue', AMI_threshold, ...
        'SMIThresholdType', 'absolute', 'SMIThresholdValue', SMI_threshold);
end


%% 8. Combined Interactive Plot: Spectrogram, AMI and SMI with Interactive Threshold Controls
disp('Creating interactive visualization...');

% Call the interactive visualization function with aperiodic data
plotInteractiveEEGVisualization(eeg_data, spect2, AMI, SMI, trimmed_time, fs, win_length, fmax, cmin, cmax, threshold_type, AMI_threshold, SMI_threshold, baseline_start_idx, baseline_end_idx, aperiodic_data, dataset_name);

disp('Interactive visualization launched. Use the controls to adjust parameters.');

%% 9. Generate Neural Derecruitment figure 

%[derecruitment_episodes, ami_derivative, time_vector] = detectDerecruitment(AMI, SMI, fs);
