function all_eeg_data = flattenConcatEDF(folderPath)
    % Default to current folder if none provided
    if nargin < 1 || isempty(folderPath)
        folderPath = '.';
    end

    % Find all .edf files in the folder
    edfFiles = dir(fullfile(folderPath, '*.edf'));
    if isempty(edfFiles)
        warning('No EDF files found in folder: %s', folderPath);
        all_eeg_data = [];
        return;
    end

    % Initialize final output array
    all_eeg_data = [];

    % Loop through each EDF file
    for i = 1:length(edfFiles)
        % Full path to the file
        filePath = fullfile(folderPath, edfFiles(i).name);

        % Read the EDF file as a timetable (hdr)
        [hdr, ~] = edfread(filePath);

        % Convert timetable to a cell array (N rows x M columns)
        % Each cell is a (samplesPerBlock x 1) vector
        eegCell = table2cell(hdr); 
        [numRows, numCols] = size(eegCell);

        % Number of samples per row-block (e.g., 1158)
        samplesPerBlock = size(eegCell{1}, 1);

        % Allocate space for flattened data from this file
        flattenedData = zeros(numRows * samplesPerBlock, numCols);

        % Flatten each column (channel) into one continuous column
        for col = 1:numCols
            % cell2mat on this column → (samplesPerBlock x numRows)
            tmp2D = cell2mat(eegCell(:, col)'); 
            % Flatten into (samplesPerBlock * numRows) x 1
            flattenedData(:, col) = tmp2D(:);
        end

        % Concatenate vertically (stack in time)
        all_eeg_data = [all_eeg_data; flattenedData];
     
    end

    % Final status
    fprintf('Done!')
end
