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
%   Last modified 03/04/2025 by Sebastian Gallo
%   galloseb@mit.edu
%
%************************************************************************** 
%% 1. Setup Environment
clear; close all; clc;
disp('Initializing State-Space Multitaper Spectrogram Analysis...');


% Select dataset (change filename to use different datasets)
datasets = {'Mass_13_Sedline_copy_raw.mat', 'SED10.mat', 'ACP_Concatenated.mat'};  % Cell array
dataset_name = datasets{3};  % Use curly braces {} to access elements

data_path = fullfile('data', dataset_name);

% Load EEG data 
data = load_eeg_data(data_path);

% Load EEG data
disp(['Loading EEG data from ', dataset_name, '...']);

fs = 178;  % Sampling frequency (Hz)
elec = 2; % Select Electrode

% Extract EEG data
if strcmp(dataset_name, 'Mass_13_Sedline_copy_raw.mat')
    eeg_data = data(1:883236, elec);
    disp(['For Black Swan #2 has been the cleanest. We only pick the ' ...
        'beginning because the rest is noise.'])
else
    % eeg_data = data(:, elec);
    eeg_data = data(1132282:2265762, elec);
end


disp('EEG data loaded successfully.');


%% 2. Preprocessing
disp('Preprocessing EEG data...');
Nt = length(eeg_data);
win_length = 2; % length of window (second)
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
%% 3. EM Algorithm
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

%% 5. Plot Spectrograms
disp('Plotting results...');

Individual_Spectrogram = true; %make false if want all spectrograms/periodograms

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
  
plot_spectrograms(spectrograms, N, win_length, fs, Individual_Spectrogram);

%% Function: Plot Spectrograms
function plot_spectrograms(spectrograms, N, win, fs, Individual_Spectrogram)
    fmax = 25;  
    cmin = -15;
    cmax = 10;

    num_plots = length(spectrograms);  % Dynamically determine number of spectrograms
    
    if Individual_Spectrogram
        fig_width = 1200;
        fig_height = 400;
        figure('Position', [100 100 fig_width fig_height]);
    else
        figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 0.7 0.9]);
    end


    for i = 1:num_plots
        subplot(num_plots, 1, i);  % Dynamically create subplots based on number of spectrograms
        title_text = spectrograms{i}{1};  % Extract subplot title
        spect_data = spectrograms{i}{2};  % Extract spectrogram data

        imagesc((1:N)*win/60, (0:fmax*win)/win, pow2db(spect_data(1:fmax*win,:)));
        axis xy;
        set(gca, 'clim', [cmin cmax]);
        ylabel('Frequency (Hz)');
        colorbar;
        colormap jet;
        title(title_text);
    end
    
    xlabel('Time (min)');
    disp('Plotting complete.');
end

