%% Parameters
fs = 125; % Sampling frequency (Hz)
window_duration = 2; % Duration of each window (seconds)
interval_duration = 50; % Duration of each interval (seconds)

%interval names
post_induction = 1;
maintenance = 2;
pre_emergence = 3;
post_extubation = 4;

% Store interval names and start times in a structured array
intervals(1) = struct('name', 'ABP 100 + CE_230', 'start_time', 230);
intervals(2) = struct('name', 'ABP 100 + CE_730', 'start_time', 730);
intervals(3) = struct('name', 'ABP 100 + CE_2600', 'start_time', 2600);
intervals(4) = struct('name', 'ABP 90 + CE_6440', 'start_time', 6440);

% Define two different start times for the two intervals (in seconds)
start_time_1 = intervals.start_time(1); % Start time for first interval
start_time_2 = 4865; % Start time for second interval (adjust as needed)

% Convert start times to sample indices
start_window_idx_1 = round(start_time_1 / window_duration);
start_window_idx_2 = round(start_time_2 / window_duration);
windows_per_interval = interval_duration / window_duration;

% Define end indices (ensuring within bounds)
end_window_idx_1 = min(start_window_idx_1 + windows_per_interval - 1, size(stateEstimate, 3));
end_window_idx_2 = min(start_window_idx_2 + windows_per_interval - 1, size(stateEstimate, 3));

% Extract state estimates for each interval
interval_stateEstimate_1 = stateEstimate(:, :, start_window_idx_1:end_window_idx_1);
interval_stateEstimate_2 = stateEstimate(:, :, start_window_idx_2:end_window_idx_2);

%% Frequency Vector
nw = size(interval_stateEstimate_1, 1); % Number of frequency bins
f = (0:(nw/2)) * (fs / nw); % Correct single-sided frequency vector
freq_mask = f <= 40; % Limit to 40 Hz
limited_f = f(freq_mask); % Apply frequency mask


%% Bootstrap Parameters
num_bootstrap = 1000; % Number of bootstrap samples
num_freq_bins = sum(freq_mask);
num_windows_1 = size(interval_stateEstimate_1, 3);
num_windows_2 = size(interval_stateEstimate_2, 3);
num_tapers = size(interval_stateEstimate_1, 2);

% Initialize storage for bootstrapped power spectra
bootstrap_power_spectra_1 = zeros(num_freq_bins, num_bootstrap);
bootstrap_power_spectra_2 = zeros(num_freq_bins, num_bootstrap);

%% Bootstrapping the Power Spectra for Each Interval
for b = 1:num_bootstrap
    % Resample windows with replacement
    resampled_indices_1 = randi(num_windows_1, [1, num_windows_1]);
    resampled_indices_2 = randi(num_windows_2, [1, num_windows_2]);

    % Storage for resampled power spectra
    resampled_power_1 = zeros(num_freq_bins, num_windows_1);
    resampled_power_2 = zeros(num_freq_bins, num_windows_2);

    % Compute bootstrapped power spectra for Interval 1
    for win_idx = 1:num_windows_1
        win = resampled_indices_1(win_idx);
        taper_power_spectra_window = zeros(num_freq_bins, num_tapers);
        for k = 1:num_tapers
            power_spectrum_k_win = abs(interval_stateEstimate_1(1:(nw/2 + 1), k, win)).^2;
            taper_power_spectra_window(:, k) = power_spectrum_k_win(freq_mask);
        end
        resampled_power_1(:, win_idx) = mean(taper_power_spectra_window, 2);
    end
    bootstrap_power_spectra_1(:, b) = pow2db(mean(resampled_power_1, 2));

    % Compute bootstrapped power spectra for Interval 2
    for win_idx = 1:num_windows_2
        win = resampled_indices_2(win_idx);
        taper_power_spectra_window = zeros(num_freq_bins, num_tapers);
        for k = 1:num_tapers
            power_spectrum_k_win = abs(interval_stateEstimate_2(1:(nw/2 + 1), k, win)).^2;
            taper_power_spectra_window(:, k) = power_spectrum_k_win(freq_mask);
        end
        resampled_power_2(:, win_idx) = mean(taper_power_spectra_window, 2);
    end
    bootstrap_power_spectra_2(:, b) = pow2db(mean(resampled_power_2, 2));
end

%% Compute Confidence Intervals
ci_lower_1 = prctile(bootstrap_power_spectra_1, 2.5, 2);
ci_upper_1 = prctile(bootstrap_power_spectra_1, 97.5, 2);
bootstrap_mean_1 = mean(bootstrap_power_spectra_1, 2);

ci_lower_2 = prctile(bootstrap_power_spectra_2, 2.5, 2);
ci_upper_2 = prctile(bootstrap_power_spectra_2, 97.5, 2);
bootstrap_mean_2 = mean(bootstrap_power_spectra_2, 2);

%% Compute Difference Between Bootstrapped Spectra

limited_f = limited_f(:); % Ensure column vector
ci_upper_1 = ci_upper_1(:); % Ensure column vector
ci_lower_1 = ci_lower_1(:); % Ensure column vector

ci_upper_2 = ci_upper_2(:); % Ensure column vector
ci_lower_2 = ci_lower_2(:); % Ensure column vector

figure;
hold on;

% Ensure vectors are the same length
fill([limited_f; flipud(limited_f)], [ci_upper_1; flipud(ci_lower_1)], 'r',...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(limited_f, bootstrap_mean_1, 'r', 'LineWidth', 1.5);

fill([limited_f; flipud(limited_f)], [ci_upper_2; flipud(ci_lower_2)], 'b',...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(limited_f, bootstrap_mean_2, 'b', 'LineWidth', 1.5);

xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Bootstrapped Power Spectra for Two Intervals');
legend('Interval 1', 'Interval 2');
grid on;
hold off;
%%
bootstrap_diff = bootstrap_power_spectra_1 - bootstrap_power_spectra_2; % Compute bootstrap differences
ci_lower_diff = prctile(bootstrap_diff, 2.5, 2); % 2.5th percentile (lower bound)
ci_upper_diff = prctile(bootstrap_diff, 97.5, 2); % 97.5th percentile (upper bound)
bootstrap_diff_mean = mean(bootstrap_diff, 2); % Mean difference


%%

ci_upper_diff = ci_upper_diff(:); % Ensure column vector
ci_lower_diff = ci_lower_diff(:); % Ensure column vector

figure;
hold on;

fill([limited_f; flipud(limited_f)], [ci_upper_diff; flipud(ci_lower_diff)],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel('Frequency (Hz)');
ylabel('Power Difference (dB)');
title('Bootstrapped Difference Between Intervals with Confidence Intervals');
grid on;
hold off;

