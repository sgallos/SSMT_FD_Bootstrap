function [spect, results] = SS_MT_cov_v3(yy, fs, TW, K, stateNoise, observationNoise, ...
                                         initialState, initialVariance)

[nw, N] = size(yy);

% Generate Slepian tapers
[tapers, ~] = dpss(nw, TW, K);

% Expand data for multitapering
y_ex = repmat(yy, 1, 1, K);
y_ex = permute(y_ex, [1, 3, 2]);

% Apply Slepian tapers
mtY = zeros(nw, K, N);
for i = 1:N
    mtY(:, :, i) = tapers .* y_ex(:, :, i);
end
mtFrequencyY = fft(mtY, nw, 1);

% Initialize Kalman filtering variables
mtStateEstimate = zeros(nw, K, N);
mtStatePrediction = zeros(nw, K, N);
mtVarianceEstimate = zeros(nw, K, N);
mtVariancePrediction = zeros(nw, K, N);
mtKalmanGain = zeros(nw, K, N);

% Full smoothed covariance matrix
mtSmoothCovariance = zeros(nw, K, K, N);

% Set initial conditions
mtStateEstimate(:, :, 1) = initialState;
mtVarianceEstimate(:, :, 1) = initialVariance;
mtSystemNoise = stateNoise;
mtObservationNoise = observationNoise;

% Forward pass: Kalman filter
for i = 2:N
    % Prediction step
    mtStatePrediction(:, :, i) = mtStateEstimate(:, :, i-1);
    mtVariancePrediction(:, :, i) = mtVarianceEstimate(:, :, i-1) + mtSystemNoise;

    % Update step
    mtKalmanGain(:, :, i) = mtVariancePrediction(:, :, i) ./ ...
        (mtObservationNoise + mtVariancePrediction(:, :, i));
    mtStateEstimate(:, :, i) = mtStatePrediction(:, :, i) + ...
        mtKalmanGain(:, :, i) .* (mtFrequencyY(:, :, i) - mtStatePrediction(:, :, i));
    mtVarianceEstimate(:, :, i) = (1 - mtKalmanGain(:, :, i)) .* mtVariancePrediction(:, :, i);
end

% Backward pass: Fixed-interval smoothing
mtSmoothState = zeros(size(mtStateEstimate));
mtSmoothVariance = zeros(size(mtVarianceEstimate));

% Initialize with forward estimates
mtSmoothState(:, :, N) = mtStateEstimate(:, :, N);
mtSmoothVariance(:, :, N) = mtVarianceEstimate(:, :, N);

for i = N-1:-1:1
    % Smoothing gain (Eq. 16)
    A_k = mtVarianceEstimate(:, :, i) ./ (mtVariancePrediction(:, :, i+1) + 1e-10);

    % Smooth state and variance (Eq. 16)
    mtSmoothState(:, :, i) = mtStateEstimate(:, :, i) + ...
        A_k .* (mtSmoothState(:, :, i+1) - mtStatePrediction(:, :, i+1));
    mtSmoothVariance(:, :, i) = mtVarianceEstimate(:, :, i) + ...
        A_k .* (mtSmoothVariance(:, :, i+1) - mtVariancePrediction(:, :, i+1)) .* A_k;

    % Full covariance matrix (Eq. 17)
    for u = 1:K
        for v = 1:K
            mtSmoothCovariance(:, u, v, i) = ...
                A_k(:, u) .* mtSmoothVariance(:, v, i+1);
        end
    end
end

% Compute multitaper spectrogram (average across tapers)
mtSpect = abs(mtSmoothState).^2;
spect = squeeze(mean(mtSpect, 2)) / fs;

% Prepare output results
results = struct(...
    'spect', spect, ...
    'stateEstimate', mtStateEstimate, ...
    'statePrediction', mtStatePrediction, ...
    'varianceEstimate', mtVarianceEstimate, ...
    'variancePrediction', mtVariancePrediction, ...
    'kalmanGain', mtKalmanGain, ...
    'smoothState', mtSmoothState, ...
    'smoothVariance', mtSmoothVariance, ...
    'smoothCovariance', mtSmoothCovariance, ...
    'systemNoise', mtSystemNoise, ...
    'observationNoise', mtObservationNoise, ...
    'mtSpect', mtSpect);

end
