% Step 1: Estimate the spectral power S_hat (⋅) using the multitaper spectral estimation of (4)
S_hat = spect2; % Ensure spect2 is computed correctly

% Step 2: Extract Fourier Coefficients
X_R_m = real(spect2_fc); % Extracts X_R^(m)(ω_j) for j = 0 to L-1
X_I_m = imag(spect2_fc); % Extracts X_I^(m)(ω_j) for j = 0 to L-1


% Get dimensions
[num_frequencies, num_tapers, num_windows] = size(X_R_m);  % Expecting (250 × 5 × 3597)

% Expand normalization_factor to match the required shape
normalization_factor = sqrt(S_hat / 2); % (250 × 3597)


% Preallocate residual matrices
s_j = zeros(num_frequencies, num_tapers, num_windows);
s_L_j = zeros(num_frequencies, num_tapers, num_windows);

% 🔹 Loop through each taper for correct element-wise division
for m = 1:num_tapers
    % Extract current taper (ensures correct shape: 250 × 3597)
    X_R_m_squeezed = squeeze(X_R_m(:, m, :));
    X_I_m_squeezed = squeeze(X_I_m(:, m, :));

    % Normalize using correctly shaped normalization_factor
    s_j(:, m, :) = X_R_m_squeezed ./ normalization_factor;
    s_L_j(:, m, :) = X_I_m_squeezed ./ normalization_factor;
end

% Concatenate along the first dimension
s_m = cat(1, s_j, s_L_j);



%% test
% Extremely important as bootstrap resample relies on assumption that
% fourier coefficients follow a normal distribution
% Compute mean and standard deviation across all elements of s_m
mean_s_m = mean(s_m(:));
std_s_m = std(s_m(:));

disp(['Mean of s_m: ', num2str(mean_s_m)])
disp(['Standard deviation of s_m: ', num2str(std_s_m)])

%% Step 3: Standardize the residuals

% Get dimensions
[num_frequencies, num_tapers, num_windows] = size(s_m); % Expecting (500 × 5 × 3597)
L = num_frequencies / 2; % Since we stored real + imaginary separately
num_total = 2 * L; % Total number of residual components

% Compute the mean of the residuals across all components for each taper & time window
mean_s_m = mean(s_m, 1); % (1 × 5 × 3597), mean over frequencies

% Subtract the mean to ensure zero mean
s_m_zero_mean = s_m - repmat(mean_s_m, [num_frequencies, 1, 1]);

% Compute standard deviation across all frequencies for each taper & time window
variance_s_m = mean(s_m_zero_mean .^ 2, 1); % Variance (1 × 5 × 3597)
std_s_m = sqrt(variance_s_m); % Standard deviation (1 × 5 × 3597)

% Normalize by the standard deviation to ensure unit variance
s_m_standardized = s_m_zero_mean ./ repmat(std_s_m, [num_frequencies, 1, 1]);

% Store the final standardized residuals
s_m = s_m_standardized;

%% Display confirmation
disp(['Mean of standardized s_m: ', num2str(mean(s_m(:)))])
disp(['Standard deviation of standardized s_m: ', num2str(std(s_m(:)))])


%% Step 4: Draw independent bootstrap replicates s_l^(m)*

num_bootstrap_samples = 1000;  % Number of bootstrap replicates
[num_frequencies, num_tapers, num_windows] = size(s_m); % Expecting (2L × M × N)
L = num_frequencies / 2; % Since real and imaginary parts are stored separately

% Initialize storage for bootstrapped spectral densities
bootstrapped_spectra = zeros(L, num_windows, num_bootstrap_samples);

% Iterate over bootstrap samples
for b = 1:num_bootstrap_samples
    % 🔹 Step 4: Draw bootstrap samples (resampling residuals with replacement)
    resampled_indices = randi([1 num_frequencies], num_frequencies, num_tapers, num_windows);

    % Resample s_m correctly using proper indexing
    s_m_bootstrap = zeros(num_frequencies, num_tapers, num_windows);
    for taper = 1:num_tapers
        for window = 1:num_windows
            s_m_bootstrap(:, taper, window) = s_m(resampled_indices(:, taper, window), taper, window);
        end
    end

    % 🔹 Step 5: Compute bootstrapped Fourier coefficients
    % Extract real and imaginary parts separately
    s_j_star = s_m_bootstrap(1:L, :, :);  % First L values: real part
    s_Lj_star = s_m_bootstrap(L+1:end, :, :);  % Last L values: imaginary part

    % Expand S_hat to match dimensions for multiplication (L × 1 × N)
    S_hat_expanded = reshape(S_hat(1:L, :), [L, 1, num_windows]);

    % Compute bootstrapped Fourier coefficients
    X_m_star = s_j_star .* sqrt(S_hat_expanded / 2) + 1i * s_Lj_star .* sqrt(S_hat_expanded / 2);

    % 🔹 Step 6: Compute bootstrapped spectral density
    S_hat_star = squeeze(mean(abs(X_m_star).^2, 2));  % Average over tapers (dimension 2)

    % Store the computed spectral density for this bootstrap iteration
    bootstrapped_spectra(:, :, b) = S_hat_star;
end

% Display confirmation
disp(['Bootstrap resampling completed: ', num2str(num_bootstrap_samples), ' replicates generated.']);


%%
% Define time range (50s–100s)
start_time = 250; % seconds
interval_duration = 50; 
end_time = start_time + interval_duration;  % seconds
window_duration = 2;

% Compute the corresponding indices in spect2
idx_start = ceil(start_time / win); % First index corresponding to >= 50s
idx_end = floor(end_time / win);    % Last index corresponding to <= 100s


% Extract bootstrapped spectra for the specified time interval
bootstrapped_spectra_trimmed = bootstrapped_spectra(:, idx_start:idx_end, :); % (250 × time_windows × 1000)

% Compute mean power spectrum across time bins
bootstrapped_spectra_mean = mean(bootstrapped_spectra_trimmed, 2); % Average over time bins
bootstrapped_spectra_mean = squeeze(bootstrapped_spectra_mean); % Remove singleton dimension

% Define confidence level (95%)
alpha = 0.05;  % Confidence level is (1 - alpha) * 100%

% Compute confidence intervals along the 3rd dimension (bootstrap samples)
ci_lower = prctile(bootstrapped_spectra_mean, 100 * (alpha / 2), 2); % 2.5th percentile
ci_upper = prctile(bootstrapped_spectra_mean, 100 * (1 - alpha / 2), 2); % 97.5th percentile

% Compute mean of the bootstrap distribution (for reference)
bootstrap_mean = mean(bootstrapped_spectra_mean, 2);


% Convert to decibels (dB)
bootstrap_mean_db = pow2db(bootstrap_mean);
ci_lower_db = pow2db(ci_lower);
ci_upper_db = pow2db(ci_upper);

% Define frequency vector
% Define correct frequency vector based on taper_power size
sf = 1 / win; % Step size in frequency
f = (0:size(bootstrapped_spectra_trimmed,1)-1) * sf; % Now matches the number of frequency bins (250)


% Plot the confidence intervals and mean power spectrum
figure;
hold on;

% Plot shaded confidence intervals
fill([f'; flipud(f')], [ci_upper_db; flipud(ci_lower_db)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot the averaged power spectrum
plot(f, bootstrap_mean_db, 'r', 'LineWidth', 1.5);

xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('Bootstrapped Power Spectrum with 95%% CI (%.2f - %.2f sec)', start_time, end_time));
grid on;
xlim([0 20]); % Adjust as needed
ylim([-10 20]); % Adjust as needed
legend('95% Confidence Interval', 'Bootstrapped Power Spectrum');
hold off;

%%
% Define time range for Interval 1
start_time_1 = 4620; % seconds
interval_duration = 50;
end_time_1 = start_time_1 + interval_duration; % seconds
window_duration = 2;

% Compute the corresponding indices for Interval 1
idx_start_1 = ceil(start_time_1 / win); % First index
idx_end_1 = floor(end_time_1 / win);   % Last index

% Extract bootstrapped spectra for Interval 1
bootstrapped_spectra_trimmed_1 = bootstrapped_spectra(:, idx_start_1:idx_end_1, :);

% Compute mean power spectrum across time bins for Interval 1
bootstrapped_spectra_mean_1 = mean(bootstrapped_spectra_trimmed_1, 2);
bootstrapped_spectra_mean_1 = squeeze(bootstrapped_spectra_mean_1);

% Define time range for Interval 2
start_time_2 = 4830; % seconds
end_time_2 = start_time_2 + interval_duration; % seconds

% Compute the corresponding indices for Interval 2
idx_start_2 = ceil(start_time_2 / win);
idx_end_2 = floor(end_time_2 / win);

% Extract bootstrapped spectra for Interval 2
bootstrapped_spectra_trimmed_2 = bootstrapped_spectra(:, idx_start_2:idx_end_2, :);

% Compute mean power spectrum across time bins for Interval 2
bootstrapped_spectra_mean_2 = mean(bootstrapped_spectra_trimmed_2, 2);
bootstrapped_spectra_mean_2 = squeeze(bootstrapped_spectra_mean_2);

% Define confidence level (95%)
alpha = 0.05; 

% Compute confidence intervals for Interval 1
ci_lower_1 = prctile(bootstrapped_spectra_mean_1, 100 * (alpha / 2), 2);
ci_upper_1 = prctile(bootstrapped_spectra_mean_1, 100 * (1 - alpha / 2), 2);

% Compute mean of the bootstrap distribution for Interval 1
bootstrap_mean_1 = mean(bootstrapped_spectra_mean_1, 2);

% Compute confidence intervals for Interval 2
ci_lower_2 = prctile(bootstrapped_spectra_mean_2, 100 * (alpha / 2), 2);
ci_upper_2 = prctile(bootstrapped_spectra_mean_2, 100 * (1 - alpha / 2), 2);

% Compute mean of the bootstrap distribution for Interval 2
bootstrap_mean_2 = mean(bootstrapped_spectra_mean_2, 2);

% Convert to decibels (dB) for both intervals
bootstrap_mean_1_db = pow2db(bootstrap_mean_1);
ci_lower_1_db = pow2db(ci_lower_1);
ci_upper_1_db = pow2db(ci_upper_1);

bootstrap_mean_2_db = pow2db(bootstrap_mean_2);
ci_lower_2_db = pow2db(ci_lower_2);
ci_upper_2_db = pow2db(ci_upper_2);

% Define frequency vector
sf = 1 / win;
f = (0:size(bootstrapped_spectra_trimmed_1,1)-1) * sf; % Matches frequency bins

% Plot the confidence intervals and mean power spectra for both intervals
figure; hold on;

% Plot shaded confidence intervals for Interval 1 (Red)
fill([f'; flipud(f')], [ci_upper_1_db; flipud(ci_lower_1_db)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot the averaged power spectrum for Interval 1 (Red)
h1 = plot(f, bootstrap_mean_1_db, 'r', 'LineWidth', 1.5);

% Plot shaded confidence intervals for Interval 2 (Blue)
fill([f'; flipud(f')], [ci_upper_2_db; flipud(ci_lower_2_db)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot the averaged power spectrum for Interval 2 (Blue)
h2 = plot(f, bootstrap_mean_2_db, 'b', 'LineWidth', 1.5);

xlabel('Frequency (Hz)', 'FontSize', 40);
ylabel('Power (dB)', 'FontSize', 40);
% title('Bootstrapped Power Spectra for Two Intervals', 'FontSize', 20);
legend([h1, h2], {'Pre-Emergence', 'Post-Extubation'}, 'FontSize', 40);

% Increase tick label size
set(gca, 'FontSize', 40);

grid on;
xlim([0 20]); % Adjust as needed
ylim([-5 30]); % Adjusted y-axis range
hold off;


