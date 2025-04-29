function [spectral_exponent, aperiodic_spectra, freqs] = extractAperiodicComponent(eeg_data, fs, win_length, overlap_pct)
% EXTRACTAPERIODICCOMPONENT - Extract aperiodic component from EEG data
%
% This function calculates the aperiodic component (1/f) from EEG data
% by fitting the power spectrum with a linear model in log-log space.
% It follows the methodology described in Widmann et al. (2025, British Journal of Anaesthesia)
%
% Inputs:
%   eeg_data      - Raw EEG data vector
%   fs            - Sampling frequency in Hz
%   win_length    - Window length in seconds (default: 2)
%   overlap_pct   - Overlap percentage between windows (default: 50)
%
% Outputs:
%   spectral_exponent - Vector of spectral exponents (alpha) over time
%   aperiodic_spectra - Matrix of aperiodic spectra (freq x windows)
%   freqs             - Frequency vector in Hz

% Set defaults if not provided
if nargin < 3
    win_length = 2; % 2-second window by default
end
if nargin < 4
    overlap_pct = 50; % 50% overlap by default
end

% Calculate window parameters in samples
window_samples = round(win_length * fs);
step_samples = round(window_samples * (1 - overlap_pct/100));
n_windows = floor((length(eeg_data) - window_samples) / step_samples) + 1;

% FFT parameters
nfft = max(2^nextpow2(window_samples), 256);

% Define frequency range (1-40 Hz)
freq_min = 1;
freq_max = 40;

% Create frequency vector within our range of interest
freq_res = fs/nfft;
freq_idx_min = ceil(freq_min/freq_res) + 1; % +1 because MATLAB is 1-indexed
freq_idx_max = floor(freq_max/freq_res) + 1;
freq_indices = freq_idx_min:freq_idx_max;
freqs = (freq_indices - 1) * freq_res; % Convert back to Hz

n_freqs = length(freqs);

% Initialize outputs
spectral_exponent = zeros(n_windows, 1);
aperiodic_spectra = zeros(n_freqs, n_windows);

% Create Hanning window for better spectral estimation
hann_window = hann(window_samples);

% Process each time window
for win_idx = 1:n_windows
    % Extract window
    start_idx = (win_idx-1) * step_samples + 1;
    end_idx = min(start_idx + window_samples - 1, length(eeg_data));
    
    % Skip if window is too small
    if end_idx - start_idx + 1 < window_samples/2
        spectral_exponent(win_idx) = NaN;
        continue;
    end
    
    % Window the data
    window_data = eeg_data(start_idx:end_idx);
    if length(window_data) < length(hann_window)
        % Pad with zeros if needed
        window_data = [window_data; zeros(length(hann_window)-length(window_data), 1)];
    end
    windowed_data = window_data .* hann_window;
    
    % Calculate power spectrum
    [pxx, f] = pwelch(windowed_data, [], [], nfft, fs);
    
    % Extract just the frequency range we want
    pxx_range = pxx(freq_indices);
    
    % Convert to log-log scale (handle zeros and negative values)
    log_freqs = log10(freqs);
    log_pxx = log10(max(pxx_range, eps)); % Avoid log(0)
    
    % Fit line to log-log spectrum
    try
        % Use polyfit instead of robustfit for simplicity and reliability
        p = polyfit(log_freqs, log_pxx, 1);
        
        % Extract slope (exponent) and offset
        exponent = p(1);
        offset = p(2);
        
        % Store spectral exponent (as positive value as in the paper)
        spectral_exponent(win_idx) = -exponent; % Convert slope to alpha
        
        % Generate aperiodic spectrum (1/f model)
        model_log_pxx = offset + exponent * log_freqs;
        aperiodic_fit = 10.^model_log_pxx;
        
        % Store in output matrix
        aperiodic_spectra(:, win_idx) = aperiodic_fit;
    catch
        % If fitting fails, mark as NaN
        spectral_exponent(win_idx) = NaN;
        aperiodic_spectra(:, win_idx) = NaN;
        warning('Fitting failed for window %d', win_idx);
    end
end

% Remove any NaN windows that might have occurred
valid_windows = ~isnan(spectral_exponent);
spectral_exponent = spectral_exponent(valid_windows);
aperiodic_spectra = aperiodic_spectra(:, valid_windows);

% Ensure we have data
if isempty(spectral_exponent)
    warning('No valid spectral exponents could be calculated. Check input data.');
    spectral_exponent = [];
    aperiodic_spectra = [];
end

end