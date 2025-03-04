%%
% Sampling frequency and interval parameters
% fs = 125;
fs = 178;
window_duration = 2; % Duration of each window (seconds)
start_time = 150*60; % Start time in seconds
interval_duration = 50; % Duration of interval (seconds)
end_time = start_time + interval_duration;

% Convert start time and duration to sample indices
start_idx_raw = round(start_time * fs);
end_idx_raw = min(start_idx_raw + interval_duration * fs - 1, length(eeg_data));

% Extract the segment of raw data for the entire interval of interest
raw_segment = eeg_data(start_idx_raw:end_idx_raw);

% Define the number of samples per window
samples_per_window_raw = window_duration * fs;

% Calculate the number of windows within the selected interval
num_windows_raw = floor(length(raw_segment) / samples_per_window_raw);

% Initialize a matrix to store power spectra for each window
P1_all_windows = zeros(floor(samples_per_window_raw / 2) + 1, num_windows_raw);

% Frequency vector for plotting
f_raw = fs * (0:(samples_per_window_raw / 2)) / samples_per_window_raw;

% Loop through each window and calculate the power spectrum
for win_raw = 1:num_windows_raw
    % Extract current window of data
    window_data = raw_segment((win_raw-1) * samples_per_window_raw + (1:samples_per_window_raw));
    
    % Compute FFT for the current window
    Y_raw = fft(window_data, samples_per_window_raw);

    % Calculate power spectrum for this window
    P2_raw = abs(Y_raw / samples_per_window_raw).^2; % Two-sided power spectrum
    P1_raw = P2_raw(1:floor(samples_per_window_raw / 2) + 1); % Single-sided power spectrum
    P1_raw(2:end-1) = 2 * P1_raw(2:end-1); % Account for negative frequencies

    % Store the power spectrum for this window
    P1_all_windows(:, win_raw) = P1_raw;
end

% Calculate the average power spectrum in the linear domain
avg_power_spectrum_raw = mean(P1_all_windows, 2);

% Convert averaged power spectrum to dB
avg_power_spectrum_db_raw = pow2db(avg_power_spectrum_raw);

% Plot individual window spectra
figure;
hold on;
for win_raw = 1:num_windows_raw
    plot(f_raw, pow2db(P1_all_windows(:, win_raw)), 'LineWidth', 0.5); % Plot each window
end
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectra of Individual Windows');
grid on;
xlim([0 20]);
ylim([-10 20]); % Adjust as needed

hold off;

% Plot averaged power spectrum
figure;
plot(f_raw, avg_power_spectrum_db_raw, 'r', 'LineWidth', 1.5); % Plot averaged spectrum
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('Average Power Spectra (%.2f - %.2f sec)', start_time, end_time))
grid on;
xlim([0 20]);
ylim([-10 20]); % Adjust as needed

%% Phenomenal validation for what is going on in a certain set of windows

% Plot each window's spectrum
figure;
hold on;
for win_raw = 1:num_windows_raw
    plot(f_raw, pow2db(P1_all_windows(:, win_raw)), 'LineWidth', 0.5); % Individual windows
end
plot(f_raw, pow2db(avg_power_spectrum_raw), 'k', 'LineWidth', 2); % Averaged spectrum
hold off;
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectra of Individual Windows and Averaged Spectrum');
grid on;
xlim([0 20]);
legend('Windows', 'Average');


