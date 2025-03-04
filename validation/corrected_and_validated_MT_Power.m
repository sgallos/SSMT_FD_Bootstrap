%%
% Sampling frequency (Hz)
% fs = 125;
fs = 178;

% Define time range (50s–100s)
start_time = 4830; % seconds
interval_duration = 50; 
end_time = start_time + interval_duration;  % seconds
window_duration = 2;

% Compute the corresponding indices in spect2
idx_start = ceil(start_time / win); % First index corresponding to >= 50s
idx_end = floor(end_time / win);    % Last index corresponding to <= 100s

idx_start = ceil(start_time / window_duration);
idx_end = floor(end_time / window_duration);

% Ensure indices are within bounds
idx_start = max(1, idx_start);
idx_end = min(size(spect2, 2), idx_end);

% Extract the portion of spect2 corresponding to 50s–100s
spect2_trimmed = spect2(:, idx_start:idx_end);
% spect2_trimmed = bootstrapped_spectra(:, idx_start:idx_end,1);

% Extract the portion of spect_taper for 50s–100s
spect_taper_trimmed = spect2_taper(:, :, idx_start:idx_end); % Frequency x Tapers x Time Windows

% Define correct frequency vector based on taper_power size
sf = 1 / win; % Step size in frequency
% sf = 1; % Step size in frequency
f = (0:size(spect_taper_trimmed,1)-1) * sf; % Now matches the number of frequency bins (250)

% Debugging check for mismatched sizes
fprintf('Corrected Frequency vector length: %d\n', length(f));
fprintf('Corrected Taper Power Size: [%d]\n', size(spect_taper_trimmed,1));

% 🎯 **Plot Individual Power Spectra for Each Taper**
figure;
hold on;
num_tapers = size(spect_taper_trimmed, 2); % Number of tapers

for taper_idx = 1:num_tapers
    % Extract power for a single taper (averaged over time bins)
    taper_power = squeeze(mean(spect_taper_trimmed(:, taper_idx, :), 3));

    % Debugging check for mismatched sizes
    fprintf('Taper %d: Length of taper_power: %d\n', taper_idx, length(taper_power));

    % Ensure taper_power has the same length as f
    if length(taper_power) ~= length(f)
        error('Mismatch between frequency vector length (%d) and taper power length (%d)', length(f), length(taper_power));
    end

    % Convert to dB and plot
    plot(f, pow2db(taper_power), 'LineWidth', 0.8);
end
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('Power Spectra of Individual Tapers (%.2f - %.2f sec)', start_time, end_time));
grid on;
xlim([0 20]); % Adjust as needed
% ylim([-10 20]); % Adjust as needed
hold off;

% 🎯 **Compute and Plot the Averaged Power Spectrum Across Tapers**
avg_power_spectrum = mean(spect2_trimmed, 2); % Average over time bins
avg_power_spectrum_db = pow2db(avg_power_spectrum); % Convert to dB

figure;
plot(f, avg_power_spectrum_db, 'r', 'LineWidth', 1.5); % Averaged power spectrum
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('Averaged Power MT Spectrum (%.2f - %.2f sec)', start_time, end_time));
grid on;
xlim([0 20]); % Adjust as needed
% ylim([-10 20]); % Adjust as needed




%%

% % Compute the corresponding indices in spect2
% idx_start = ceil(start_time / win); % First index corresponding to >= 30s
% idx_end = floor(end_time / win);    % Last index corresponding to <= 80s
% 
% % Ensure indices are within bounds
% idx_start = max(1, idx_start);
% idx_end = min(size(spect2, 2), idx_end);
% 
% % Extract the portion of spect2 corresponding to 30s–80s
% spect2_trimmed = spect2(:, idx_start:idx_end);
% 
% % Define the new time vector (in minutes, to match imagesc)
% time_trimmed = (idx_start:idx_end) * win;
% 
% % Define the frequency vector (Hz)
% f = (0:fmax*win) * sf;
% 
% % Plot the extracted spectrogram
% figure;
% imagesc(time_trimmed, f, pow2db(spect2_trimmed(1:fmax*win, :)));
% axis xy;
% set(gca, 'clim', [cmin cmax]);
% ylabel('Frequency (Hz)');
% xlabel('Time (sec)');
% colorbar;
% title(sprintf('Multitaper Spectrogram (%.2f - %.2f sec)', start_time, end_time));
% drawnow;

