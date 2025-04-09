%% Simulate AR(6) process (same as before)
% Sampling rate and signal length
Fs = 178;
duration_sec = 100;
N = Fs * duration_sec;

% AR(6) coefficients from Kim et al.
ar_coeffs = [3.9515, -7.8885, 9.7340, -7.7435, 3.8078, -0.9472];

% Simulate AR signal
v = randn(N,1);
x = zeros(N,1);
for t = 7:N
    x(t) = ar_coeffs * flipud(x(t-6:t-1)) + v(t);
end


%% Define Chronux Parameters and Compute Chronux Multitaper Spectrogram

% Define tapers: TW = 3, K = 3
params.Fs = 178;
params.tapers = [3, 3];
params.pad = 0;  % no zero-padding
params.fpass = [0, params.Fs/2];  % full Nyquist range
params.err = [0 0];
params.trialave = 0;  % do not average over trials

movingwin = [2, 0.5];  % [window size, step size] in seconds

[S_chronux, t_chronux, f_chronux] = mtspecgramc(x, movingwin, params);
S_chronux_dB = pow2db(S_chronux');  % Transpose for freq × time

% Compute color limits for contrast
clim = prctile(S_chronux_dB(:), [5 95]);

figure;
imagesc(t_chronux/60, f_chronux, S_chronux_dB, clim);  % time in minutes
axis xy;
colormap jet;

xlabel('Time (Minutes)', 'FontSize', 16);
ylabel('Frequency (Hz)', 'FontSize', 16);
title('Chronux Multitaper Spectrogram', 'FontSize', 18);
c = colorbar('southoutside');
ylabel(c, 'Power (dB)', 'FontSize', 14);
set(gca, 'FontSize', 14);


%% Compute Manual Multitaper Spectrogram
Fs = 178;
N = length(x);              % total signal length
win_len_sec = 2;            % window length in seconds
step_sec = 0.5;             % step size in seconds
TW = 3;                     % time-bandwidth product
K = 3;                      % number of tapers

win_len = round(win_len_sec * Fs);
step = round(step_sec * Fs);

% Time vector (center of each window)
t_centers = (1:step:(N - win_len + 1)) + floor(win_len/2);
t_sec = t_centers / Fs;

% Get tapers once
[tapers, ~] = dpss(win_len, TW, K);

% Preallocate output: freq × time
nFreqs = win_len/2 + 1;
nTime = length(t_sec);
S_dB = zeros(nFreqs, nTime);

for i = 1:nTime
    idx = (i-1)*step + (1:win_len);
    x_win = x(idx);

    J = zeros(nFreqs, K);  % FFTs per taper
    for k = 1:K
        tapered = x_win .* tapers(:,k);
        fft_vals = fft(tapered);
        J(:,k) = fft_vals(1:nFreqs);
    end

    % Compute average power
    S = mean(abs(J).^2, 2);
    S_dB(:, i) = pow2db(S);
end

% Frequency axis
f_Hz = linspace(0, Fs/2, nFreqs);

% Compute color limits from percentiles for good contrast
clim = prctile(S_dB(:), [5 95]);

figure;
imagesc(t_sec/60, f_Hz, S_dB, clim);  % time in minutes (or just t_sec)
axis xy;
colormap jet;

xlabel('Time (Minutes)', 'FontSize', 16);
ylabel('Frequency (Hz)', 'FontSize', 16);
title('Manual Multitaper Spectrogram', 'FontSize', 18);
c = colorbar('southoutside');
ylabel(c, 'Power (dB)', 'FontSize', 14);
set(gca, 'FontSize', 14);
%%

Fs = 178;
N = 1024;  % or 5 min × Fs for longer signal

% AR(6) coefficients from paper
ar_coeffs = [3.9515, -7.8885, 9.7340, -7.7435, 3.8078, -0.9472];
a = [1, -ar_coeffs];

% Generate AR(6) signal with unit-variance white noise
x = filter(1, a, randn(N,1));  % Do NOT normalize

% DPSS tapers
TW = 3;
K = 3;
[tapers, ~] = dpss(N, TW, K);

% FFT of tapered signals
J = zeros(N/2+1, K);  % positive frequencies only
for k = 1:K
    x_tapered = x .* tapers(:,k);
    Xk = fft(x_tapered);
    J(:,k) = Xk(1:N/2+1);  % one-sided
end

% Multitaper power spectral density
S_mt = mean(abs(J).^2, 2);           % power
S_mt_dB = pow2db(S_mt);              % convert to dB
f_mt = linspace(0, Fs/2, N/2+1);     % frequency axis
f_mt_norm = f_mt / Fs;               % normalized freq (0 to 0.5)


[H, f_true] = freqz(1, a, N/2+1, Fs);  % match resolution
S_true = abs(H).^2;
S_true_dB = pow2db(S_true);
f_true_norm = f_true / Fs;


%% Perform Monte Carlo "Gold Standard", Part 1: Plot the true spectrum alongside the manual and chronux estimate
[S_chronux, f_chronux] = mtspectrumc(x, params);
S_chronux_dB = pow2db(S_chronux);
f_chronux_norm = f_chronux / Fs;

figure; hold on;

% True spectrum
plot(f_true_norm, S_true_dB, 'g-', 'LineWidth', 2);

% Manual multitaper
plot(f_mt_norm, S_mt_dB, 'k', 'LineWidth', 2);

% Chronux multitaper
plot(f_chronux_norm, S_chronux_dB, 'r--', 'LineWidth', 2);

% Labels, title, legend
xlabel('Normalized Frequency (Hz / Fs)', 'FontSize', 16);
ylabel('Power (dB)', 'FontSize', 16);
legend('True Spectrum', 'Manual Multitaper', 'Chronux Multitaper', ...
       'Location', 'northeast');
title('True vs. Manual vs. Chronux Multitaper Spectrum (AR(6))');
xlim([0 0.5]);
grid on;
set(gca, 'FontSize', 14);

%% Perform Monte Carlo "Gold Standard", Part 2: (Paper-Faithful: Full-Length Manual Multitaper)
% Based on Kim et al., Fig. 1 — 1000 trials of AR(6), full-length MT spectrum

% Parameters (from paper)
nTrials = 1000;
N = 1024;           % Full signal length used in the paper
Fs = 178;           % Sampling rate
TW = 3; K = 3;      % Time-bandwidth and taper count

% Frequency axis
nFreqs = N/2 + 1;
f = linspace(0, Fs/2, nFreqs);
f_norm = f / Fs;

%% EQUATION 1 FROM KIM et. al 2018 PAPER
% DPSS tapers (computed once)
[tapers, ~] = dpss(N, TW, K);
%%

% Preallocate power matrix (in linear scale)
spectra_lin = zeros(nTrials, nFreqs);

% Loop over Monte Carlo simulations
for i = 1:nTrials
    % Simulate AR(6) signal (white noise input with unit variance)
    x_i = filter(1, [1, -ar_coeffs], randn(N,1));

    % Manual multitaper spectrum (single window, K tapers)
    J = zeros(nFreqs, K);
    for k = 1:K
        x_tapered = x_i .* tapers(:,k);
        fft_vals = fft(x_tapered);
        J(:,k) = fft_vals(1:nFreqs);  % keep positive frequencies
    end

    % Average across tapers → store linear power
    S_i = mean(abs(J).^2, 2);
    spectra_lin(i,:) = S_i;

    if mod(i, 100) == 0
        fprintf('Completed %d / %d simulations...\n', i, nTrials);
    end
end

% Save raw results
save('MC_manual_MT_spectra_linear.mat', 'spectra_lin', 'f', 'Fs');

% Compute Confidence Intervals and Plot
% Compute empirical confidence intervals in **linear power**
CI_lower_lin = prctile(spectra_lin, 2.5);
CI_upper_lin = prctile(spectra_lin, 97.5);
MT_mean_lin  = mean(spectra_lin, 1);

% Convert to dB after aggregation (faithful to paper’s approach)
CI_lower_dB = pow2db(CI_lower_lin);
CI_upper_dB = pow2db(CI_upper_lin);
MT_mean_dB  = pow2db(MT_mean_lin);

%% Select One Manual Multitaper Estimate (Single Realization)
% Simulate one AR(6) realization
x_single = filter(1, [1, -ar_coeffs], randn(N,1));

% Compute multitaper spectrum manually
J_single = zeros(nFreqs, K);
for k = 1:K
    x_tapered = x_single .* tapers(:,k);
    fft_vals = fft(x_tapered);
    J_single(:,k) = fft_vals(1:nFreqs);
end
S_single = mean(abs(J_single).^2, 2);
S_single_dB = pow2db(S_single);

% Plot: MC CI + True + Single MT Estimate
figure; hold on;

% Confidence interval (shaded region)
fill([f_norm, fliplr(f_norm)], ...
     [CI_lower_dB, fliplr(CI_upper_dB)], ...
     [0.85 0.9 1], 'EdgeColor', 'none');

% True spectrum (green)
plot(f_true_norm, S_true_dB, 'g-', 'LineWidth', 2);

% Single multitaper estimate (black)
plot(f_norm, S_single_dB, 'k', 'LineWidth', 2);

xlabel('Normalized Frequency (Hz / Fs)', 'FontSize', 16);
ylabel('Power (dB)', 'FontSize', 16);
legend('95% CI (MC)', 'True Spectrum', 'One MT Estimate', 'Location', 'northeast');
title('Monte Carlo Confidence Interval with One MT Estimate');
xlim([0 0.5]);
grid on;
set(gca, 'FontSize', 14);

%% Jackknife procedure

% Parameters (from your pipeline)
N = 1024;
Fs = 178;
TW = 3;
K = 3;
nFreqs = N/2 + 1;
f = linspace(0, Fs/2, nFreqs);
f_norm = f / Fs;

% DPSS tapers
[tapers, ~] = dpss(N, TW, K);

% Simulate one AR(6) time series
x = filter(1, [1, -ar_coeffs], randn(N,1));  % white noise input, unit variance


%% EQUATION 3 FROM KIM et. al 2018 PAPER
% Compute tapered FFTs (J_k)
J = zeros(nFreqs, K);
for k = 1:K
    x_tapered = x .* tapers(:,k);
    Xk = fft(x_tapered);
    J(:,k) = Xk(1:nFreqs);  % keep positive freqs
end
%% EQUATION 2 FROM KIM et. al 2018 PAPER


% Individual taper spectra (linear power)
S_k = abs(J).^2;  % size: [nFreqs × K]

%% Preallocate: one jackknife replicate per taper
S_jack = zeros(K, nFreqs);  % [K × nFreqs]

for k = 1:K
    idx = setdiff(1:K, k);  % indices of tapers to include
    S_jack(k, :) = mean(S_k(:, idx), 2);  % leave-one-out average
end

%% Step 1: Log-transform the jackknife spectra
log_S_jack = log(S_jack);  % natural log (not dB)

% Step 2: Compute jackknife mean log-spectrum
log_S_jack_mean = mean(log_S_jack, 1);  % size: [1 × nFreqs]

% Step 3: Jackknife variance estimate at each frequency
JK_var = (K - 1)/K * sum((log_S_jack - log_S_jack_mean).^2, 1);  % [1 × nFreqs]

% Step 4: Compute critical t-value
alpha = 0.05;
t_crit = tinv(1 - alpha/2, K - 1);  % two-sided 95% CI, df = K-1

% Step 5: Compute CI bounds (in log space)
log_CI_upper = log_S_jack_mean + t_crit * sqrt(JK_var);
log_CI_lower = log_S_jack_mean - t_crit * sqrt(JK_var);

% Step 6: Convert to dB
CI_lower_dB_JK = 10 * log10(exp(log_CI_lower));
CI_upper_dB_JK = 10 * log10(exp(log_CI_upper));

%% EQUATION 4 FROM KIM et. al 2018 PAPER: Compute multitaper estimate from full K tapers (same signal)
S_mt = mean(S_k, 2);           % full estimate in linear power
%%

S_mt_dB = pow2db(S_mt);        % convert to dB

% Plot
figure; hold on;x

% Monte Carlo 95% CI (shaded region)
fill([f_norm, fliplr(f_norm)], ...
     [CI_lower_dB, fliplr(CI_upper_dB)], ...
     [0.85 0.9 1], 'EdgeColor', 'none');  % pale blue fill

% True spectrum
plot(f_norm, S_true_dB, 'g-', 'LineWidth', 2);

% Manual MT estimate (black)
plot(f_norm, S_mt_dB, 'k-', 'LineWidth', 2);

% Jackknife confidence interval bounds
plot(f_norm, CI_lower_dB_JK, 'm--', 'LineWidth', 1.5);  % magenta = lower
plot(f_norm, CI_upper_dB_JK, 'c--', 'LineWidth', 1.5);  % cyan = upper

% Labels and formatting
xlabel('Normalized Frequency (Hz / Fs)', 'FontSize', 16);
ylabel('Power (dB)', 'FontSize', 16);
legend('95% CI (Monte Carlo)', ...
       'True Spectrum', ...
       'Multitaper Estimate', ...
       'Jackknife Lower', ...
       'Jackknife Upper', ...
       'Location', 'northeast');
title('Jackknife and Monte Carlo Confidence Intervals (Manual MT)', 'FontSize', 18);
xlim([0 0.5]);
grid on;
set(gca, 'FontSize', 14);

%% MFDB

% --- J is the [nFreqs × K] matrix of tapered Fourier transforms ---
% It was computed from:
% for k = 1:K
%     x_tapered = x .* tapers(:,k);
%     fft_vals = fft(x_tapered);
%     J(:,k) = fft_vals(1:nFreqs);
% end

% STEP 2: MFDB Bootstrapping from J_k(f)

nBoot = 1000;  % number of bootstrap replicates
S_mfdb = zeros(nBoot, nFreqs);  % [nBoot × freqs] for storing power spectra

%% EQUATIONS 7 & 8 FROM KIM et. al 2018 PAPER
for b = 1:nBoot
    % Sample with replacement from K tapers
    resample_idx = randi(K, [1, K]);

    % Resample columns of J based on those indices
    J_resampled = J(:, resample_idx);

    % Compute multitaper spectrum (mean power across resampled J's)
    S_boot = mean(abs(J_resampled).^2, 2);  % linear power
    S_mfdb(b, :) = S_boot;
end
%%

% Compute 95% confidence intervals in linear space
CI_lower_mfdb_lin = prctile(S_mfdb, 2.5);
CI_upper_mfdb_lin = prctile(S_mfdb, 97.5);

% Convert to dB
CI_lower_dB_mfdb = pow2db(CI_lower_mfdb_lin);
CI_upper_dB_mfdb = pow2db(CI_upper_mfdb_lin);

figure; hold on;

% True spectrum (green)
plot(f_norm, S_true_dB, 'g-', 'LineWidth', 2);

% Multitaper estimate (black)
plot(f_norm, S_mt_dB, 'k-', 'LineWidth', 2);

% Jackknife confidence intervals (magenta / cyan dashed)
plot(f_norm, CI_lower_dB_JK, 'c--', 'LineWidth', 1.5);  % lower
plot(f_norm, CI_upper_dB_JK, 'm--', 'LineWidth', 1.5);  % upper

% MFDB confidence intervals (red / blue dashed)
plot(f_norm, CI_lower_dB_mfdb, 'r--', 'LineWidth', 2.5);  % lower
plot(f_norm, CI_upper_dB_mfdb, 'b--', 'LineWidth', 2.5);  % upper

% Formatting
xlabel('Normalized Frequency (Hz / Fs)', 'FontSize', 16);
ylabel('Power (dB)', 'FontSize', 16);
legend('True Spectrum', ...
       'Multitaper Estimate', ...
       'Jackknife Lower', 'Jackknife Upper', ...
       'MFDB Lower', 'MFDB Upper', ...
       'Location', 'northeast');
title('MFDB vs. Jackknife Confidence Intervals (Manual MT)', 'FontSize', 18);
xlim([0 0.5]);
grid on;
set(gca, 'FontSize', 14);

