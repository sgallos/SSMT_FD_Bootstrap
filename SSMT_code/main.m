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

close all; clear all; clc;

disp('Loading data...');
subjName = 'SED10'; 
subj = [subjName,'.mat']; 
% data = load(subj);

% data=data.data;
% fs = 250; %Sampling rate BU 
% fs = 178.5; %Sampling rate Black Swan
fs = 178; %Sampling rate validation
% fs = 125; %Sampling rate Swine
channel=1; % Frontal Channel
% 
% Get meaningful data from EEG data
% y = data(channel,600*fs+1:1880*fs)';
% load('Mass_13_Sedline_copy_filtered.mat');
% load('Mass_13_Sedline_Lowpass_20.mat');
load('Mass_13_Sedline_copy_raw.mat');
% % eeg_data = data_Mass_13_Sedline_copy_raw(:,2);
eeg_data = data_Mass_13_Sedline_copy_raw(1:883236,2);



% load('EEG_data_run_2 copy.mat'); %%%%%USING FOR SWINE EXPERIMENTS
% eeg_data = EEG_waveform;
% EEG_time_adjusted = EEG_time - 1.7316e9 + 41.9840;
% idx_start = find(EEG_time_adjusted >= 2890, 1, 'first');
% eeg_data = EEG_waveform(1, idx_start:end);

% load('validation_F8.mat');
% eeg_data = electrodeData;
% eeg_data = electrodeData(520000:2615000,:);
% load('EEG_waveform_run_3.mat')
% eeg_data = data_Mass_13_Sedline_copy_filtered(1:882000,2);

% %%
% fs = 178; % Sampling rate (Hz)
% [b, a] = butter(4, [.5 20]/(fs/2), 'bandpass'); % 4th order Butterworth filter
% eeg_data = filtfilt(b, a, eeg_data);


%%

% load('SED10.mat');
% eeg_data= data(1,:);


% eeg_data = EEG_waveform;
% y = eeg_data(600*fs+1:1880*fs);
y = eeg_data;
% y = data(1,:);

Nt = length(y);
% clear data;

win = 2; % length of window (second)
sf = 1/win; % one step of frequency
nw = win*fs; % the number of elements in a window
N = floor(Nt/nw); % the number of window

% Matrix form of data according to the size of window
yy = reshape(y(1:nw*N),nw,N); 

rtaper = rectwin(nw);
rtaper = rtaper/sqrt(nw);

kyy = zeros(nw,N);
for i = 1 : N
    kyy(:,i) = rtaper.*yy(:,i);
end
% Fourier transfrom of data
frequencyY = fft(kyy,nw,1);

%% EM Algorithm
% First we can limit the frequency range to 0 to 30 Hz and we can adjust
% the max level depnding on the EEG data for greater denoising. 
OBSNOISE_CUTOFF = 30*win; % 30 Hz

% Initially we can set the alpha and beta as 1 
alpha = 1;
beta = 1;

% Initial guess for the observation noise
observationNoise = 100;
GUESS_WINDOW_LENGTH = 150; % EM estimatin: 5 min = 300 sec = 300/win = 150
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

%% Spectral estimation (periodogram, multitaper, SS-P, SS-MT)

% periodogram
spect1 = periodogram(yy, fs);
% multitaper spectrogram
[spect2,spect2_taper, spect2_fc] = multitaper_fc(yy, fs, TW, K);
% state-space periodogram
spect3 = SS_ST(yy, fs, sn, on, is, iv);
% statate-space multitaper spectrogram
% [spect4, results_MT] = SS_MT(yy, fs, TW, K, mtSn, mtOn, mtIs, mtIv);
[spect4, results_MT] = SS_MT_cov_v3(yy, fs, TW, K, mtSn, mtOn, mtIs, mtIv);
% [spect4, results_MT] = SS_MT_cov_v2(yy, fs, TW, K, mtSn, mtOn, mtIs, mtIv);

%% Plot for spectrogram comparisons

fig = figure('color','w','units','normalized','position',[0 0 0.7 0.9]); clf;

fmax = 25;
cmin = -15;
cmax = 10;

colormap jet

ax(1) = subplot(411);
imagesc((1:N)*win/60, (0:fmax*win)*sf, pow2db(spect1(1:fmax*win,:)));
axis xy;
set(gca,'clim',[cmin cmax])
set(ax(1),'xticklabel',[]);
ylabel('Frequency (Hz)');
colorbar
title('Periodogram')
drawnow

ax(2) = subplot(412);
imagesc((1:N)*win/60, (0:fmax*win)*sf, pow2db(spect2(1:fmax*win,:)));
axis xy;
set(gca,'clim',[cmin cmax])
set(ax(2),'xticklabel',[]);
ylabel('Frequency (Hz)');
colorbar
title('Multitaper Spectrogram')
drawnow

ax(3) = subplot(413);
imagesc( (1:N)*win/60, (0:fmax*win)*sf, pow2db(spect3(1:fmax*win,:)));
axis xy;
set(gca,'clim',[cmin cmax])
set(ax(3),'xticklabel',[]);
ylabel('Frequency (Hz)');
colorbar
title('State-Space Periodogram')
drawnow

ax(4) = subplot(414);
imagesc( (1:N)*win/60, (0:fmax*win)*sf, pow2db(spect4(1:fmax*win,:)));
axis xy;
set(gca,'clim',[cmin cmax])
ylabel('Frequency (Hz)');
xlabel('Time (min)');
colorbar
title('State-Space Multitaper Spectrogram')
drawnow

%%
% % Adjust EEG time to remove the large offset
% EEG_time_adjusted = EEG_time - 1.7316e9 + 41.9840;
% 
% % Find the index where EEG_time_adjusted reaches 2824.28
% idx_start = find(EEG_time_adjusted >= 2824.28, 1, 'first');
% 
% % Trim EEG time to start at 2824.28
% EEG_time_trimmed = EEG_time_adjusted(idx_start:end);
% spect4_trimmed = spect4(:, idx_start:end);  % Trim spectrogram data
% 
% % Ensure t_spect is correctly spaced and aligned with EEG_time_trimmed
% t_spect_4 = EEG_time_trimmed;  % Use the actual EEG time values

%%
figure;
imagesc(((1:N)*win)/60, (0:fmax*win)*sf, pow2db(spect4(1:fmax*win,:)));
axis xy;
colormap jet;
set(gca, 'clim', [-15 10]); % Color limits
set(gca, 'FontSize', 15);
xlabel('Time (sec)', 'FontSize', 15);
ylabel('Frequency (Hz)', 'FontSize', 15);
h = colorbar; % Add color bar
ylabel(h, 'Power (dB)', 'FontSize', 15); % Label the color bar
% title('State-Space Multitaper Spectrogram');

%%
% spect2_trimmed = spect2(:, idx_start:end);  % Trim spectrogram data
% 
% % Ensure t_spect is correctly spaced and aligned with EEG_time_trimmed
% t_spect_2 = EEG_time_trimmed;  % Use the actual EEG time values


%%
% figure;
% % imagesc(((1:N)*win)/60, (0:fmax*win)*sf, pow2db(spect4(1:fmax*win,:)));
% imagesc(t_spect_2, (0:fmax*win)*sf, pow2db(spect2(1:fmax*win,:)));
% axis xy;
% colormap jet;
% set(gca, 'clim', [-15 10]); % Color limits
% xlabel('Time (sec)');
% ylabel('Frequency (Hz)');
% colorbar;
% title('State-Space Multitaper Spectrogram');

%%
figure;
imagesc((1:N)*win, (0:fmax*win)*sf, pow2db(spect2(1:fmax*win,:)));
axis xy;
colormap jet;
set(gca, 'clim', [-15 10]); % Color limits
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
h = colorbar; % Add color bar
ylabel(h, 'Power (dB)'); % Label the color bar
title('State-Space Multitaper Spectrogram');

%% Code for obtaining the posterior mean & posterior variance for the joint posterior distribution
% Sebastian Gallo, 11/04/2024

stateEstimate = results_MT.stateEstimate;
varianceEstimate = results_MT.varianceEstimate;
% covarianceMatrix = results_MT.covarianceMatrix;
covarianceMatrix = results_MT.smoothCovariance;

% %% Plot EEG
% % Define the sampling frequency
% % fs = 178.5; % Sampling frequency in Hz
% fs = 178; % Sampling frequency in Hz
% 
% % Generate a time vector based on the EEG data length
% time = (0:length(eeg_data)-1) / fs; % Time in seconds
% 
% % Plot the filtered EEG data
% figure;
% plot(time, eeg_data);
% xlabel('Time (seconds)', 'FontSize', 12);
% ylabel('Amplitude (\muV)', 'FontSize', 12);
% title('Filtered EEG Data', 'FontSize', 14);
% grid on;
