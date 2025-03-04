function eeg_data = load_eeg_data(filename)
    % Function to load all EEG data from a .mat file
    % It automatically detects the variable containing 'data' in its name.
    %
    % Usage:
    %   eeg_data = load_eeg_data('data/SED10.mat');  % Loads all EEG data
    
    % Load the .mat file into a structure
    data = load(filename);
    
    % Get all variable names inside the .mat file
    fieldNames = fieldnames(data);
    
    % Look for a variable that contains "data" in its name
    matches = contains(fieldNames, 'data', 'IgnoreCase', true);
    
    % If one or more matches are found, pick the first one
    if any(matches)
        selectedVar = fieldNames{find(matches, 1)}; % Select first match
    else
        warning('No variable with "data" found in %s. Using first available variable: %s', filename, fieldNames{1});
        selectedVar = fieldNames{1}; % Fallback to first variable
    end

    % Extract the entire EEG matrix
    eeg_data = data.(selectedVar);

    % Display confirmation message
    fprintf('Loaded EEG data from %s (Variable: %s, Size: %dx%d)\n', ...
        filename, selectedVar, size(eeg_data, 1), size(eeg_data, 2));
end

