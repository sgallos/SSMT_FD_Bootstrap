function downsampleRippleData(inputFile, newSamplingRate, outputFile)
% DOWNSAMPLERIPPLEDATA imports data from Ripple Neuromed files, downsamples it, and saves as a .mat file
%
% Inputs:
%   inputFile - Full path to the input file (.nf3 or .nev or .ns*)
%   newSamplingRate - Target sampling rate in Hz (e.g., 250)
%   outputFile - Full path for the output .mat file
%
% Need to add "neuroshare" folder and subfolder to path in order to handle
% the nf3 data.
%
% Example usage:
%   downsampleRippleData('C:\Users\chris\Documents\EEG Analysis\Data\datafile.nf3', 250, 'C:\Users\chris\Documents\EEG Analysis\Data\downsampleddata.mat');

fprintf('Starting EEG data downsampling process...\n');

% Open file using Neuroshare API
fprintf('Opening file %s...\n', inputFile);
[ns_RESULT, hFile] = ns_OpenFile(inputFile);

if ~strcmp(ns_RESULT, 'ns_OK')
    error('Error opening file: %s', ns_RESULT);
end

% Get file information
[ns_RESULT, fileInfo] = ns_GetFileInfo(hFile);
if ~strcmp(ns_RESULT, 'ns_OK')
    ns_CloseFile(hFile);
    error('Error getting file info: %s', ns_RESULT);
end

% Find all analog entities (continuous data channels)
analogEntities = [];
channelLabels = {};
originalSamplingRate = 0;

fprintf('Finding analog data channels...\n');
for i = 1:fileInfo.EntityCount
    [~, entityInfo] = ns_GetEntityInfo(hFile, i);
    
    % Check if entity is analog (continuous) data
    if entityInfo.EntityType == 2 % Analog entity
        analogEntities = [analogEntities, i];
        channelLabels{end+1} = entityInfo.EntityLabel;
        
        % Get the sampling rate (get it from the first analog channel)
        if originalSamplingRate == 0
            [~, analogInfo] = ns_GetAnalogInfo(hFile, i);
            originalSamplingRate = analogInfo.SampleRate;
        end
    end
end

if isempty(analogEntities)
    ns_CloseFile(hFile);
    error('No analog data channels found in the file');
end

fprintf('Found %d analog channels at %d Hz\n', length(analogEntities), originalSamplingRate);

% Calculate downsampling parameters
downsampleFactor = round(originalSamplingRate / newSamplingRate);
actualNewRate = originalSamplingRate / downsampleFactor;

fprintf('Downsampling from %d Hz to %.2f Hz (factor: %d)\n', ...
    originalSamplingRate, actualNewRate, downsampleFactor);

% Get the number of data points
[~, entityInfo] = ns_GetEntityInfo(hFile, analogEntities(1));
numSamples = entityInfo.ItemCount;

% Process data in chunks to avoid memory issues
chunkSizeSec = 60; % Process 60 seconds at a time
chunkSizeSamples = min(chunkSizeSec * originalSamplingRate, numSamples);
numChunks = ceil(numSamples / chunkSizeSamples);

% Initialize output matrix - as TIME x CHANNELS (samples as rows, channels as columns)
outputLength = ceil(numSamples / downsampleFactor);
downsampledData = zeros(outputLength, length(analogEntities));

% Set up anti-aliasing filter
nyquist = originalSamplingRate / 2;
cutoff = 0.8 * (newSamplingRate / 2); % 80% of new Nyquist frequency
[b, a] = butter(8, cutoff / nyquist, 'low');

% Process data in chunks
fprintf('Processing data in %d chunks...\n', numChunks);
tic

for chunk = 1:numChunks
    fprintf('Processing chunk %d of %d (%.1f%% complete)...\n', ...
        chunk, numChunks, 100 * (chunk-1) / numChunks);
    
    % Calculate start and end indices for this chunk
    startIdx = (chunk - 1) * chunkSizeSamples + 1;
    endIdx = min(chunk * chunkSizeSamples, numSamples);
    chunkSize = endIdx - startIdx + 1;
    
    % Read data for all channels in this chunk
    for ch = 1:length(analogEntities)
        entityID = analogEntities(ch);
        
        % Get analog data for this chunk
        [ns_RESULT, ~, chunkData] = ns_GetAnalogData(hFile, entityID, startIdx, chunkSize);
        
        if ~strcmp(ns_RESULT, 'ns_OK')
            warning('Error reading data for channel %d: %s', ch, ns_RESULT);
            chunkData = zeros(chunkSize, 1);
        end
        
        % Apply anti-aliasing filter
        filteredData = filtfilt(b, a, double(chunkData));
        
        % Downsample
        downsampledChunk = filteredData(1:downsampleFactor:end);
        
        % Calculate output indices
        outStartIdx = ceil(startIdx / downsampleFactor);
        outEndIdx = outStartIdx + length(downsampledChunk) - 1;
        
        % Store in output array - note we're storing in column ch
        if outEndIdx <= size(downsampledData, 1)
            downsampledData(outStartIdx:outEndIdx, ch) = downsampledChunk;
        end
    end
end

% Create header information
header = struct();
header.originalFile = inputFile;
header.originalSamplingRate = originalSamplingRate;
header.downsampledRate = actualNewRate;
header.downsampleFactor = downsampleFactor;
header.channelLabels = channelLabels;
header.numChannels = length(analogEntities);
header.filterInfo = struct('type', 'butterworth', 'order', 8, 'cutoff', cutoff);

% Add dimension info to header
header.dimensions = 'samples x channels'; % Explicitly noting the data organization

% Save downsampled data
fprintf('Saving data to %s\n', outputFile);
save(outputFile, 'downsampledData', 'header', '-v7.3');

% Close the file
ns_CloseFile(hFile);

% Report time taken
processingTime = toc;
fprintf('Processing completed in %.1f seconds (%.2f minutes)\n', ...
    processingTime, processingTime/60);
fprintf('Downsampled data saved successfully!\n');
end