%% Parameters
fs = 125; % Sampling frequency (Hz)
window_duration = 2; % Duration of each window (seconds)
interval_duration = 50; % Duration of each interval (seconds)

% Store interval names and start times in a structured array
intervals(1) = struct('name', 'ABP 100 + CE_(230)', 'start_time', 230);
intervals(2) = struct('name', 'ABP 100 + CE_(730)', 'start_time', 730);
intervals(3) = struct('name', 'ABP 100 + CE_(2600)', 'start_time', 2600);
intervals(4) = struct('name', 'ABP 90 + CE_(6440)', 'start_time', 6440);

num_intervals = length(intervals);

%% Frequency Vector
nw = size(stateEstimate, 1); % Number of frequency bins
f = (0:(nw/2)) * (fs / nw); % Correct single-sided frequency vector
freq_mask = f <= 20; % Limit to 20 Hz
limited_f = f(freq_mask); % Apply frequency mask

%% Bootstrap Parameters
num_bootstrap = 1000; % Number of bootstrap samples
num_freq_bins = sum(freq_mask); % 41 Total Hz Bins
num_tapers = size(stateEstimate, 2); % 3 Tapers

%% Loop Over Unique Interval Pairs
for i = 1:num_intervals
    for j = i+1:num_intervals  % Avoid duplicate comparisons

        % Extract interval names and start times
        name_1 = intervals(i).name;
        name_2 = intervals(j).name;
        start_time_1 = intervals(i).start_time;
        start_time_2 = intervals(j).start_time;

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

        % Number of windows
        num_windows_1 = size(interval_stateEstimate_1, 3);
        num_windows_2 = size(interval_stateEstimate_2, 3);

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

        %% Plot Bootstrapped Power Spectra
        
        limited_f = limited_f(:); % Ensure column vector
        ci_upper_1 = ci_upper_1(:); % Ensure column vector
        ci_lower_1 = ci_lower_1(:); % Ensure column vector
        
        ci_upper_2 = ci_upper_2(:); % Ensure column vector
        ci_lower_2 = ci_lower_2(:); % Ensure column vector







        figure;
        hold on;
        
        % Plot confidence intervals (shaded areas)
        fill([limited_f; flipud(limited_f)], [ci_upper_1; flipud(ci_lower_1)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill([limited_f; flipud(limited_f)], [ci_upper_2; flipud(ci_lower_2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % Plot mean power spectra and store handles
        h1 = plot(limited_f, bootstrap_mean_1, 'r', 'LineWidth', 1.5);
        h2 = plot(limited_f, bootstrap_mean_2, 'b', 'LineWidth', 1.5);
        
        % Correctly assign legend to the mean plots only
        legend([h1, h2], name_1, name_2);
        
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        title([name_1 ' vs ' name_2]);
        grid on;
        hold off;

        %% Compute and Plot Bootstrapped Difference
        bootstrap_diff = bootstrap_power_spectra_1 - bootstrap_power_spectra_2;
        ci_lower_diff = prctile(bootstrap_diff, 2.5, 2);
        ci_upper_diff = prctile(bootstrap_diff, 97.5, 2);
        bootstrap_diff_mean = mean(bootstrap_diff, 2);
        ci_upper_diff = ci_upper_diff(:); % Ensure column vector
        ci_lower_diff = ci_lower_diff(:); % Ensure column vector

        figure;
        hold on;
        
        % Plot confidence interval (shaded area)
        h_fill = fill([limited_f; flipud(limited_f)], [ci_upper_diff; flipud(ci_lower_diff)], ...
                      'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % Add a solid black horizontal line at y = 0
        yline(0, 'k', 'LineWidth', 1.5);
        
        % Create an invisible line for the legend
        h_legend = plot(nan, nan, 'k', 'LineWidth', 1.5);
        
        % Assign the legend correctly
        legend(h_legend, '95% Confidence Interval');
        
        xlabel('Frequency (Hz)');
        ylabel('Power Difference (dB)');
        title(['Bootstrapped Difference: ' name_1 ' vs ' name_2]);
        grid on;
        hold off;

    end
end
