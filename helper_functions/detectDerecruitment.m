function [ami_episodes, smi_episodes, ami_derivative, smi_derivative, time_vector] = detectDerecruitment(ami, smi, fs, drop_threshold, window_minutes, do_plot)
% DETECTDERECRUITMENT - Detect neural derecruitment episodes from AMI and SMI
%
% This function detects periods of neural derecruitment by analyzing the 
% rate of change of Alpha Modulation Index (AMI) and Slow Modulation Index (SMI).
%
% Inputs:
%   ami           - Vector containing Alpha Modulation Index time series
%   smi           - Vector containing Slow Modulation Index time series
%   fs            - Sampling frequency of the indices in Hz
%   drop_threshold - Rate of decrease threshold (default: 0.1 per minute)
%   window_minutes - Time window for analysis in minutes (default: 1)
%   do_plot       - Whether to plot results (default: true)
%
% Outputs:
%   ami_episodes  - Nx2 matrix with start and end indices of AMI derecruitment episodes
%   smi_episodes  - Nx2 matrix with start and end indices of SMI derecruitment episodes
%   ami_derivative - Smoothed AMI derivative in units of AMI change per minute
%   smi_derivative - Smoothed SMI derivative in units of SMI change per minute
%   time_vector   - Time vector in minutes

% Set default parameters if not provided
if nargin < 4
    drop_threshold = 0.1;
end
if nargin < 5
    window_minutes = 1;
end
if nargin < 6
    do_plot = true;
end

% Make sure ami and smi are the same length
if length(ami) ~= length(smi)
    error('AMI and SMI vectors must be the same length');
end

% Create time vector in minutes
time_vector = (0:length(ami)-1)/(fs*60);

% Calculate window size in samples
window_samples = round(window_minutes * 60 * fs);

% Check if vectors are too short
if length(ami) < window_samples
    error('Input vectors are shorter than the analysis window. Please provide more data or use a shorter window.');
end

% Process AMI
ami_episodes = processModulationIndex(ami, fs, drop_threshold, window_minutes, 'AMI', time_vector);

% Process SMI
smi_episodes = processModulationIndex(smi, fs, drop_threshold, window_minutes, 'SMI', time_vector);

% Calculate derivatives for output
ami_derivative = calculateDerivative(ami, fs);
smi_derivative = calculateDerivative(smi, fs);

% Plot results if requested
if do_plot
    plotDerecruitmentResults(ami, smi, time_vector, fs, ami_episodes, smi_episodes, drop_threshold, window_minutes);
end

end % End of main function

% Function to process modulation index and find episodes
function episodes = processModulationIndex(signal, fs, drop_threshold, window_minutes, signal_name, time_vector)
    % Calculate derivative using a centered difference formula
    derivative_raw = zeros(size(signal));
    for i = 2:length(signal)-1
        derivative_raw(i) = (signal(i+1) - signal(i-1)) / (2/fs);
    end
    
    % Handle endpoints
    derivative_raw(1) = (signal(2) - signal(1)) / (1/fs);
    derivative_raw(end) = (signal(end) - signal(end-1)) / (1/fs);
    
    % Convert derivative to change per minute
    derivative_raw = derivative_raw * 60;
    
    % Apply smoothing using a moving average to reduce noise
    window_size = round(fs * 5); % 5-second smoothing window
    if mod(window_size, 2) == 0
        window_size = window_size + 1; % Make sure window size is odd
    end
    derivative_smooth = movmean(derivative_raw, window_size);
    
    % Detect derecruitment episodes (where derivative is below negative threshold)
    is_derecruitment = derivative_smooth < -drop_threshold;
    
    % Find the start and end of each derecruitment episode
    state_changes = diff([0; is_derecruitment; 0]);
    starts = find(state_changes == 1);
    ends = find(state_changes == -1) - 1;
    
    % Combine starts and ends
    episodes = [starts, ends];
    
    if ~isempty(episodes)
        % Calculate actual drop for each episode
        signal_drops = zeros(size(episodes, 1), 1);
        
        for i = 1:size(episodes, 1)
            start_idx = episodes(i, 1);
            end_idx = episodes(i, 2);
            signal_drops(i) = signal(start_idx) - signal(end_idx);
            
            % Calculate duration in minutes
            duration_mins = (end_idx - start_idx) / (fs * 60);
            
            % Adjust end index for episodes that are shorter than window_minutes
            if duration_mins < window_minutes
                potential_end_idx = min(start_idx + round(window_minutes * 60 * fs), length(signal));
                if potential_end_idx > end_idx
                    potential_drop = signal(start_idx) - signal(potential_end_idx);
                    % If we still get enough drop by extending to the full window, update end index
                    if potential_drop >= drop_threshold
                        end_idx = potential_end_idx;
                        signal_drops(i) = potential_drop;
                        episodes(i, 2) = end_idx;
                    end
                end
            end
        end
        
        % Only keep episodes with drop of at least threshold over the time window
        valid_episodes = signal_drops >= drop_threshold;
        episodes = episodes(valid_episodes, :);
    end
    
    % Filter out very short episodes (less than 3 seconds)
    if ~isempty(episodes)
        min_duration = round(3 * fs); % 3 seconds in samples
        valid_episodes = (episodes(:,2) - episodes(:,1) + 1) >= min_duration;
        episodes = episodes(valid_episodes, :);
    end
    
    % Print only a summary to console, not individual episodes
    fprintf('Detected %d %s derecruitment episodes\n', size(episodes, 1), signal_name);
end

% Function to calculate smoothed derivatives
function derivative = calculateDerivative(signal, fs)
    % Calculate derivative using a centered difference formula
    derivative_raw = zeros(size(signal));
    for i = 2:length(signal)-1
        derivative_raw(i) = (signal(i+1) - signal(i-1)) / (2/fs);
    end
    
    % Handle endpoints
    derivative_raw(1) = (signal(2) - signal(1)) / (1/fs);
    derivative_raw(end) = (signal(end) - signal(end-1)) / (1/fs);
    
    % Convert derivative to change per minute
    derivative_raw = derivative_raw * 60;
    
    % Apply smoothing using a moving average to reduce noise
    window_size = round(fs * 5); % 5-second smoothing window
    if mod(window_size, 2) == 0
        window_size = window_size + 1; % Make sure window size is odd
    end
    derivative = movmean(derivative_raw, window_size);
end

% Function to plot results
function plotDerecruitmentResults(ami, smi, time_vector, fs, ami_episodes, smi_episodes, drop_threshold, window_minutes)
    % Create a figure with two subplots for AMI and SMI
    fig = figure('Position', [100, 100, 1200, 700], 'Name', 'Neural Derecruitment Detection');
    
    % Store data in figure for callback access
    fig.UserData.ami = ami;
    fig.UserData.smi = smi;
    fig.UserData.time_vector = time_vector;
    fig.UserData.fs = fs;
    fig.UserData.drop_threshold = drop_threshold;
    fig.UserData.window_minutes = window_minutes;
    
    % Plot 1: AMI with derecruitment episodes highlighted
    ax1 = subplot(2, 1, 1);
    plot(time_vector, ami, 'LineWidth', 1.5);
    hold on;
    
    % Highlight AMI derecruitment episodes
    highlightEpisodes(ax1, ami, time_vector, ami_episodes, fs, 'r');
    
    title('AMI with Detected Derecruitment Episodes', 'FontSize', 12);
    ylabel('Alpha Modulation Index (AMI)');
    grid on;
    xlim([min(time_vector), max(time_vector)]);
    ylim([max(0, min(ami)-0.1), min(1, max(ami)+0.1)]);
    set(gca, 'XTickLabel', []); % Remove x-axis labels from top plot
    
    % Add legend
    legend('AMI', 'Derecruitment Episodes', 'Location', 'southeast');
    
    % Plot 2: SMI with derecruitment episodes highlighted
    ax2 = subplot(2, 1, 2);
    plot(time_vector, smi, 'LineWidth', 1.5);
    hold on;
    
    % Highlight SMI derecruitment episodes
    highlightEpisodes(ax2, smi, time_vector, smi_episodes, fs, 'r');
    
    title('SMI with Detected Derecruitment Episodes', 'FontSize', 12);
    xlabel('Time (minutes)');
    ylabel('Slow Modulation Index (SMI)');
    grid on;
    xlim([min(time_vector), max(time_vector)]);
    ylim([max(0, min(smi)-0.1), min(1, max(smi)+0.1)]);
    
    % Add legend
    legend('SMI', 'Derecruitment Episodes', 'Location', 'southeast');
    
    % Add summary text for AMI
    [ami_summary_text, ~] = generateSummaryText(ami, ami_episodes, fs, 'AMI');
    ami_summary_box = annotation('textbox', [0.75, 0.75, 0.2, 0.15], 'String', ami_summary_text, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % Add summary text for SMI
    [smi_summary_text, ~] = generateSummaryText(smi, smi_episodes, fs, 'SMI');
    smi_summary_box = annotation('textbox', [0.75, 0.35, 0.2, 0.15], 'String', smi_summary_text, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % Add overall title
    sgtitle('Neural Derecruitment Detection', 'FontWeight', 'bold', 'FontSize', 14);
    
    % Link the x axes of the two plots
    linkaxes([ax1, ax2], 'x');
    
    % Store handle references for updates
    fig.UserData.ax1 = ax1;
    fig.UserData.ax2 = ax2;
    fig.UserData.ami_summary_box = ami_summary_box;
    fig.UserData.smi_summary_box = smi_summary_box;
    
    % Add interactive controls
    % Create a panel for the controls
    control_panel = uipanel('Title', 'Detection Parameters', 'Position', [0.01, 0.01, 0.25, 0.08]);
    
    % Add text labels
    uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Index Drop:', ...
        'Position', [10, 35, 60, 20], 'HorizontalAlignment', 'left');
    uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Time Window (min):', ...
        'Position', [130, 35, 120, 20], 'HorizontalAlignment', 'left');
    
    % Add edit boxes for threshold and window
    threshold_edit = uicontrol('Parent', control_panel, 'Style', 'edit', ...
        'String', num2str(drop_threshold), 'Position', [70, 15, 50, 20], ...
        'Callback', @updateParameters);
    
    window_edit = uicontrol('Parent', control_panel, 'Style', 'edit', ...
        'String', num2str(window_minutes), 'Position', [250, 15, 50, 20], ...
        'Callback', @updateParameters);
    
    % Store the edit box handles
    fig.UserData.threshold_edit = threshold_edit;
    fig.UserData.window_edit = window_edit;
    
    % Add a button to export the results
    export_button = uicontrol('Style', 'pushbutton', 'String', 'Export Results', ...
        'Position', [350, 15, 100, 25], 'Callback', @exportResults);
    
    % Add zoom/pan buttons
    zoom_button = uicontrol('Style', 'pushbutton', 'String', 'Zoom', ...
        'Position', [460, 15, 60, 25], 'Callback', @toggleZoom);
    
    pan_button = uicontrol('Style', 'pushbutton', 'String', 'Pan', ...
        'Position', [530, 15, 60, 25], 'Callback', @togglePan);
    
    reset_view_button = uicontrol('Style', 'pushbutton', 'String', 'Reset View', ...
        'Position', [600, 15, 80, 25], 'Callback', @resetView);
    
    % Callback function to update the plot when parameters change
    function updateParameters(source, ~)
        % Get the current values from the edit boxes
        new_threshold = str2double(fig.UserData.threshold_edit.String);
        new_window = str2double(fig.UserData.window_edit.String);
        
        % Validate inputs
        if isnan(new_threshold) || new_threshold <= 0
            new_threshold = fig.UserData.drop_threshold;
            fig.UserData.threshold_edit.String = num2str(new_threshold);
            warning('Invalid threshold value. Using previous value.');
        end
        
        if isnan(new_window) || new_window <= 0
            new_window = fig.UserData.window_minutes;
            fig.UserData.window_edit.String = num2str(new_window);
            warning('Invalid window value. Using previous value.');
        end
        
        % Update stored values
        fig.UserData.drop_threshold = new_threshold;
        fig.UserData.window_minutes = new_window;
        
        % Re-detect derecruitment episodes with new parameters
        [new_ami_episodes, new_smi_episodes] = detectDerecruitment(fig.UserData.ami, fig.UserData.smi, ...
            fig.UserData.fs, new_threshold, new_window, false);
        
        % Clear previous plot content
        cla(fig.UserData.ax1);
        cla(fig.UserData.ax2);
        
        % Redraw AMI plot
        plot(fig.UserData.ax1, fig.UserData.time_vector, fig.UserData.ami, 'LineWidth', 1.5);
        hold(fig.UserData.ax1, 'on');
        
        % Redraw AMI derecruitment episodes
        highlightEpisodes(fig.UserData.ax1, fig.UserData.ami, fig.UserData.time_vector, new_ami_episodes, fig.UserData.fs, 'r');
        hold(fig.UserData.ax1, 'off');
        legend(fig.UserData.ax1, 'AMI', 'Derecruitment Episodes', 'Location', 'southeast');
        
        % Redraw SMI plot
        plot(fig.UserData.ax2, fig.UserData.time_vector, fig.UserData.smi, 'LineWidth', 1.5);
        hold(fig.UserData.ax2, 'on');
        
        % Redraw SMI derecruitment episodes
        highlightEpisodes(fig.UserData.ax2, fig.UserData.smi, fig.UserData.time_vector, new_smi_episodes, fig.UserData.fs, 'r');
        hold(fig.UserData.ax2, 'off');
        legend(fig.UserData.ax2, 'SMI', 'Derecruitment Episodes', 'Location', 'southeast');
        
        % Update summary texts
        [ami_summary_text, ~] = generateSummaryText(fig.UserData.ami, new_ami_episodes, fig.UserData.fs, 'AMI');
        fig.UserData.ami_summary_box.String = ami_summary_text;
        
        [smi_summary_text, ~] = generateSummaryText(fig.UserData.smi, new_smi_episodes, fig.UserData.fs, 'SMI');
        fig.UserData.smi_summary_box.String = smi_summary_text;
        
        % Print summary to console
        fprintf('\nUpdated parameters: Index drop = %.2f, Window = %.2f min\n', new_threshold, new_window);
    end
    
    % Function to export results
    function exportResults(~, ~)
        % Create a dialog to get filename
        [filename, pathname] = uiputfile({'*.mat', 'MAT-files (*.mat)'; ...
                                          '*.csv', 'CSV-files (*.csv)'}, ...
                                         'Save Derecruitment Results As');
        if isequal(filename, 0) || isequal(pathname, 0)
            return;
        end
        
        [~, ~, ext] = fileparts(filename);
        
        % Re-detect derecruitment episodes with current parameters to get fresh results
        [ami_episodes, smi_episodes, ami_derivative, smi_derivative] = detectDerecruitment(fig.UserData.ami, fig.UserData.smi, ...
            fig.UserData.fs, fig.UserData.drop_threshold, fig.UserData.window_minutes, false);
        
        % Generate summary data
        [~, ami_summary] = generateSummaryText(fig.UserData.ami, ami_episodes, fig.UserData.fs, 'AMI');
        [~, smi_summary] = generateSummaryText(fig.UserData.smi, smi_episodes, fig.UserData.fs, 'SMI');
        
        % Save data based on file extension
        if strcmpi(ext, '.mat')
            % Save as MAT file
            drop_threshold = fig.UserData.drop_threshold; % Variable for saving
            window_minutes = fig.UserData.window_minutes; % Variable for saving
            
            save(fullfile(pathname, filename), 'ami_episodes', 'smi_episodes', ...
                'ami_derivative', 'smi_derivative', 'ami_summary', 'smi_summary', ...
                'drop_threshold', 'window_minutes');
            fprintf('Results saved to %s\n', fullfile(pathname, filename));
        end
    end
    
    % Function to toggle zoom mode
    function toggleZoom(~, ~)
        zoom on;
        pan off;
    end
    
    % Function to toggle pan mode
    function togglePan(~, ~)
        zoom off;
        pan on;
    end
    
    % Function to reset view
    function resetView(~, ~)
        zoom off;
        pan off;
        
        % Reset axes to show all data
        xlim(fig.UserData.ax1, [min(fig.UserData.time_vector), max(fig.UserData.time_vector)]);
        ylim(fig.UserData.ax1, [max(0, min(fig.UserData.ami)-0.1), min(1, max(fig.UserData.ami)+0.1)]);
        
        ylim(fig.UserData.ax2, [max(0, min(fig.UserData.smi)-0.1), min(1, max(fig.UserData.smi)+0.1)]);
    end
end

% Helper function to highlight episodes
function highlightEpisodes(ax, signal, time, episodes, fs, color)
    for i = 1:size(episodes, 1)
        start_idx = episodes(i, 1);
        end_idx = episodes(i, 2);
        
        % Calculate time points
        episode_time = time(start_idx:end_idx);
        episode_signal = signal(start_idx:end_idx);
        
        % Plot episode
        plot(ax, episode_time, episode_signal, color, 'LineWidth', 2.5);
        
        % No text annotations as requested
    end
end

% Helper function to generate summary text
function [summary_text, summary_data] = generateSummaryText(signal, episodes, fs, signal_name)
    num_episodes = size(episodes, 1);
    summary_data = struct('num_episodes', num_episodes);
    
    if num_episodes > 0
        % Calculate total derecruitment time
        total_seconds = sum((episodes(:,2) - episodes(:,1)))/fs;
        total_minutes = total_seconds/60;
        percent_time = 100 * total_seconds / (length(signal)/fs);
        
        % Calculate average drop
        avg_drop = 0;
        for i = 1:size(episodes, 1)
            start_idx = episodes(i, 1);
            end_idx = episodes(i, 2);
            avg_drop = avg_drop + (signal(start_idx) - signal(end_idx));
        end
        avg_drop = avg_drop / num_episodes;
        
        summary_data.total_minutes = total_minutes;
        summary_data.percent_time = percent_time;
        summary_data.avg_drop = avg_drop;
        
        % Simplified summary text without total minutes and percentage
        summary_text = sprintf('%s Derecruitment:\n%d episodes\nAvg drop: %.2f', ...
            signal_name, num_episodes, avg_drop);
    else
        summary_text = sprintf('No %s derecruitment\nepisodes detected', signal_name);
        summary_data.total_minutes = 0;
        summary_data.percent_time = 0;
        summary_data.avg_drop = 0;
    end
end