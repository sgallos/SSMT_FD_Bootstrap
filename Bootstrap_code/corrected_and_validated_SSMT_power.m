% Sampling frequency (Hz)
fs = 178;

% Define time range (50s–100s)
start_time = 250; % seconds
interval_duration = 50; 
end_time = start_time + interval_duration;  % seconds

% Compute the corresponding indices in the spectrogram
idx_start = ceil(start_time / win); % First index corresponding to >= 50s
idx_end = floor(end_time / win);    % Last index corresponding to <= 100s

% Ensure indices are within bounds
idx_start = max(1, idx_start);
idx_end = min(size(results_MT.spect, 2), idx_end);

% Extract the portion of state-space spectrogram for 50s–100s
SS_spect_trimmed = results_MT.spect(:, idx_start:idx_end);
SS_spect_taper_trimmed = results_MT.mtSpect(:, :, idx_start:idx_end); % Frequency x Tapers x Time Windows

% Define frequency vector using spectrogram logic
sf = 1 / win; % Step size in frequency
f = (0:size(SS_spect_taper_trimmed,1)-1) * sf; % Now matches the number of frequency bins

% Debugging check for mismatched sizes
fprintf('Corrected Frequency vector length: %d\n', length(f));
fprintf('Corrected Taper Power Size: [%d]\n', size(SS_spect_taper_trimmed,1));

% 🎯 **Plot Individual Power Spectra for Each Taper**
figure;
hold on;
num_tapers = size(SS_spect_taper_trimmed, 2); % Number of tapers

for taper_idx = 1:num_tapers
    % Extract power for a single taper (averaged over time bins)
%     taper_power = squeeze(mean(SS_spect_taper_trimmed(:, taper_idx, :), 3));
      taper_power = squeeze(sum(SS_spect_taper_trimmed(:, taper_idx, :), 3) / size(SS_spect_taper_trimmed,3)) / fs;


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
title(sprintf('State-Space Power Spectra of Individual Tapers (%.2f - %.2f sec)', start_time, end_time));
grid on;
xlim([0 20]); % Adjust as needed
ylim([-10 20]); % Adjust as needed
hold off;

% 🎯 **Compute and Plot the Averaged Power Spectrum Across Tapers**
avg_power_spectrum = mean(SS_spect_trimmed, 2); % Average over time bins
avg_power_spectrum_db = pow2db(avg_power_spectrum); % Convert to dB

figure;
plot(f, avg_power_spectrum_db, 'r', 'LineWidth', 1.5); % Averaged power spectrum
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('Averaged State-Space Power Spectrum (%.2f - %.2f sec)', start_time, end_time));
grid on;
xlim([0 20]); % Adjust as needed
ylim([-10 20]); % Adjust as needed