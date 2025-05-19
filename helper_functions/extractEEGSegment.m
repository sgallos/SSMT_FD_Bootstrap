function extractEEGSegment(eeg_data, trimmed_time, fs, start_time, end_time, output_filename)
% EXTRACTEEGSEGMENT - Extract a specific time window from EEG data and save it
%
% This function extracts a segment of EEG data between specified time points
% and saves it in a format compatible with main_v3.m
%
% Usage:
%   extractEEGSegment(eeg_data, trimmed_time, fs, start_time, end_time, output_filename)
%
% Inputs:
%   eeg_data        - Raw EEG data vector
%   trimmed_time    - Time vector in minutes corresponding to eeg_data
%   fs              - Sampling frequency in Hz
%   start_time      - Start time in minutes
%   end_time        - End time in minutes
%   output_filename - Base filename for export (without extension)
%
% Example:
%   extractEEGSegment(eeg_data, trimmed_time, fs, 10.5, 15.2, 'segment1');

% Input validation
if nargin < 6 || isempty(output_filename)
    output_filename = 'extracted_eeg_segment';
end

% Validate time range
if start_time >= end_time
    error('Start time must be less than end time');
end

% Keep times within data bounds
start_time = max(min(trimmed_time), min(start_time, max(trimmed_time)));
end_time = max(min(trimmed_time), min(end_time, max(trimmed_time)));

% Find indices corresponding to selected time window
start_idx = find(trimmed_time >= start_time, 1, 'first');
end_idx = find(trimmed_time <= end_time, 1, 'last');

if isempty(start_idx) || isempty(end_idx) || start_idx >= end_idx
    error('Invalid time range for the provided data');
end

% Extract the EEG data within the selected time range
extracted_eeg = eeg_data(start_idx:end_idx);

% Format the data to match main_v3.m expectations - single column matrix
data = extracted_eeg(:);

% Ask user where to save the file
[file, path] = uiputfile('*.mat', 'Save extracted data as', [output_filename '.mat']);

if file ~= 0
    % Save the file with the variable named 'data' as required by main_v3.m
    save(fullfile(path, file), 'data');
    
    % Create a companion info file with metadata
    info = struct();
    info.extraction_date = datestr(now);
    info.start_time_min = start_time;
    info.end_time_min = end_time;
    info.duration_min = end_time - start_time;
    info.duration_sec = (end_time - start_time) * 60;
    info.fs = fs;
    info.samples = length(extracted_eeg);
    
    % Save the info file
    [~, base_name, ~] = fileparts(file);
    info_file = [base_name '_info.mat'];
    save(fullfile(path, info_file), 'info');
    
    % Display success message
    fprintf('Extracted %d samples (%.2f min) from %.2f to %.2f min\n', ...
        length(extracted_eeg), end_time-start_time, start_time, end_time);
    fprintf('Data saved to: %s\n', fullfile(path, file));
    fprintf('Info saved to: %s\n', fullfile(path, info_file));
end
end