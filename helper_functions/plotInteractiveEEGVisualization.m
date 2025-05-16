function plotInteractiveEEGVisualization(eeg_data, spect2, AMI, SMI, trimmed_time, fs, win_length, fmax, cmin, cmax, threshold_type, AMI_threshold, SMI_threshold, baseline_start_idx, baseline_end_idx, aperiodic_data)
% PLOTINTERACTIVEEEGVISUALIZATION Creates an interactive figure with spectrogram, AMI, and SMI plots
%   This function creates a figure with interactive controls for adjusting threshold
%   types and values, as well as spectrogram color scale parameters.
%
% Inputs:
%   eeg_data - The raw EEG data
%   spect2 - The precomputed spectrogram data
%   AMI - Precomputed Alpha Modulation Index
%   SMI - Precomputed Slow Modulation Index
%   trimmed_time - Time vector (in minutes)
%   fs - Sampling frequency
%   win_length - Window length in seconds
%   fmax - Maximum frequency to display
%   cmin - Minimum color scale value (dB)
%   cmax - Maximum color scale value (dB)
%   threshold_type - 'baseline' or 'absolute'
%   AMI_threshold - Initial AMI threshold value (for absolute mode)
%   SMI_threshold - Initial SMI threshold value (for absolute mode)
%   baseline_start_idx - Start index of baseline period (seconds)
%   baseline_end_idx - End index of baseline period (seconds)

% Create a figure with 4 subplots (added raw EEG plot)
fig = figure('Name', 'EEG Analysis', 'Position', [50, 50, 1500, 900], 'units', 'normalized', 'MenuBar', 'none', 'ToolBar', 'figure');

% Store data in the figure's UserData for access by the callbacks
fig.UserData.eeg_data = eeg_data;
fig.UserData.fs = fs;
fig.UserData.baseline_start_idx = baseline_start_idx;
fig.UserData.baseline_end_idx = baseline_end_idx;
fig.UserData.trimmed_time = trimmed_time;
fig.UserData.current_threshold_type = threshold_type;
fig.UserData.AMI_threshold = AMI_threshold;
fig.UserData.SMI_threshold = SMI_threshold;
fig.UserData.cmin = cmin;
fig.UserData.cmax = cmax;
fig.UserData.spectrogram_data = spect2(1:fmax*win_length,:);
fig.UserData.AMI = AMI;
fig.UserData.SMI = SMI;
fig.UserData.aperiodic_data = aperiodic_data;
fig.UserData.show_raw_eeg = false; % Initialize flag for raw EEG display

% Create tabs for different view modes
tabgroup = uitabgroup('Parent', fig, 'Position', [0, 0, 1, 1]);
tab1 = uitab(tabgroup, 'Title', 'Standard View');
tab2 = uitab(tabgroup, 'Title', 'EEG + Spectrogram View');
tab3 = uitab(tabgroup, 'Title', '1/f Component View');

% ===== TAB 1: STANDARD VIEW =====
% 1. Spectrogram subplot
ax1 = axes('Parent', tab1, 'Position', [0.05, 0.60, 0.75, 0.35]);
h_spec = imagesc(trimmed_time, (0:fmax*win_length)/win_length, pow2db(spect2(1:fmax*win_length,:)));
axis xy;
set(gca, 'clim', [cmin cmax]);
title('Multitaper Spectrogram', 'FontWeight', 'bold');
ylabel('Frequency (Hz)');
% Add colorbar to the right with proper positioning
cb = colorbar('Position', [0.81, 0.60, 0.015, 0.35]);
cb.Label.String = 'Power (dB)';
colormap(jet);
set(gca, 'XTickLabel', []); % Remove x-tick labels for all but the bottom subplot

% Store spectrogram handle for updating
fig.UserData.h_spec = h_spec;
fig.UserData.cb = cb;
fig.UserData.ax1 = ax1;

% 2. Plot AMI - medium height ratio
ax2 = axes('Parent', tab1, 'Position', [0.05, 0.325, 0.75, 0.225]);
ami_plot = plot(trimmed_time, AMI, 'LineWidth', 1.5);
title('Alpha Modulation Index (AMI)', 'FontWeight', 'bold');
ylabel('AMI');
ylim([0 1]);
grid on;
set(gca, 'XTickLabel', []); % Remove x-tick labels for all but the bottom subplot

% 3. SMI subplot - medium height ratio
ax3 = axes('Parent', tab1, 'Position', [0.05, 0.05, 0.75, 0.225]);
smi_plot = plot(trimmed_time, SMI, 'LineWidth', 1.5);
title('Slow Modulation Index (SMI)', 'FontWeight', 'bold');
xlabel('Time (min)');
ylabel('SMI');
ylim([0 1]);
grid on;

% Store the plot handles for updating later
fig.UserData.ami_plot = ami_plot;
fig.UserData.smi_plot = smi_plot;
fig.UserData.ax2 = ax2;
fig.UserData.ax3 = ax3;

% Link all x-axes so zooming in one affects all
linkaxes([ax1, ax2, ax3], 'x');

% ===== TAB 2: EEG + SPECTROGRAM VIEW =====
% 1. Raw EEG subplot
ax_eeg = axes('Parent', tab2, 'Position', [0.05, 0.60, 0.75, 0.35]);
% Calculate reasonable time vector for raw EEG (depending on sampling frequency)
eeg_time_min = linspace(trimmed_time(1), trimmed_time(end), length(eeg_data));
raw_eeg_plot = plot(ax_eeg, eeg_time_min, eeg_data, 'LineWidth', 0.5);
title('Raw EEG Signal', 'FontWeight', 'bold');
ylabel('Amplitude (μV)');
grid on;
set(ax_eeg, 'XTickLabel', []); % Remove x-tick labels

% 2. Spectrogram subplot beneath raw EEG
ax_spec2 = axes('Parent', tab2, 'Position', [0.05, 0.05, 0.75, 0.5]);
h_spec2 = imagesc(ax_spec2, trimmed_time, (0:fmax*win_length)/win_length, pow2db(spect2(1:fmax*win_length,:)));
axis xy;
set(ax_spec2, 'clim', [cmin cmax]);
title('Multitaper Spectrogram', 'FontWeight', 'bold');
xlabel('Time (min)');
ylabel('Frequency (Hz)');
% Add colorbar to the right
cb2 = colorbar('Position', [0.81, 0.05, 0.015, 0.5]);
cb2.Label.String = 'Power (dB)';

% Store the raw EEG handles
fig.UserData.raw_eeg_plot = raw_eeg_plot;
fig.UserData.ax_eeg = ax_eeg;
fig.UserData.h_spec2 = h_spec2;
fig.UserData.ax_spec2 = ax_spec2;
fig.UserData.cb2 = cb2;

% Link the x-axes on tab2
linkaxes([ax_eeg, ax_spec2], 'x');

% BANDPASS FILTER ADDITIONS TO EEG VISUALIZATION

% Store original EEG data for reset purposes
fig.UserData.original_eeg_data = eeg_data;
fig.UserData.original_spectrogram = spect2(1:fmax*win_length,:);

% Adjust the existing plot positions in tab2 to make room for power plot
set(fig.UserData.ax_eeg, 'Position', [0.05, 0.70, 0.75, 0.25]);   % Raw EEG at top
set(fig.UserData.ax_spec2, 'Position', [0.05, 0.40, 0.75, 0.25]); % Spectrogram in middle

% Create the power trajectory plot at the bottom
ax_power = axes('Parent', tab2, 'Position', [0.05, 0.10, 0.75, 0.25]);
power_plot = plot(ax_power, trimmed_time, zeros(size(trimmed_time)), 'LineWidth', 1.5);
title(ax_power, 'Bandpass Filtered Power', 'FontWeight', 'bold');
ylabel(ax_power, 'Power');
xlabel(ax_power, 'Time (min)');
grid(ax_power, 'on');

% Store the power plot handles
fig.UserData.ax_power = ax_power;
fig.UserData.power_plot = power_plot;

% Link all three axes
linkaxes([fig.UserData.ax_eeg, fig.UserData.ax_spec2, ax_power], 'x');

% ===== TAB 3: APERIODIC (1/f) COMPONENT VIEW =====
% 1. Original spectrogram subplot at top
ax_spec3 = axes('Parent', tab3, 'Position', [0.05, 0.70, 0.75, 0.25]);
h_spec3 = imagesc(ax_spec3, trimmed_time, (0:fmax*win_length)/win_length, pow2db(spect2(1:fmax*win_length,:)));
axis(ax_spec3, 'xy');
set(ax_spec3, 'clim', [cmin cmax]);
title(ax_spec3, 'Original Spectrogram', 'FontWeight', 'bold');
ylabel(ax_spec3, 'Frequency (Hz)');
cb3 = colorbar('Position', [0.81, 0.70, 0.015, 0.25]);
cb3.Label.String = 'Power (dB)';
set(ax_spec3, 'XTickLabel', []);

% 2. Spectral exponent (alpha) plot in middle
ax_alpha = axes('Parent', tab3, 'Position', [0.05, 0.40, 0.75, 0.25]);
h_alpha = plot(ax_alpha, aperiodic_data.exponent_time, aperiodic_data.spectral_exponent, 'LineWidth', 1.5);
title(ax_alpha, 'Spectral Exponent (\alpha)', 'FontWeight', 'bold');
ylabel(ax_alpha, '\alpha');
grid(ax_alpha, 'on');
set(ax_alpha, 'XTickLabel', []);

% Add consciousness threshold markers based on paper values
hold(ax_alpha, 'on');
awake_threshold = 1.80;  % Upper bound of awake range
anesthetized_threshold = 2.36;  % Lower bound of anesthetized range
plot(ax_alpha, [min(trimmed_time), max(trimmed_time)], [awake_threshold, awake_threshold], '--', 'Color', [0, 0.7, 0], 'LineWidth', 1.5);
plot(ax_alpha, [min(trimmed_time), max(trimmed_time)], [anesthetized_threshold, anesthetized_threshold], '--', 'Color', [0.7, 0, 0], 'LineWidth', 1.5);
text(ax_alpha, max(trimmed_time)*0.02, awake_threshold*1.02, 'Awake Threshold', 'Color', [0, 0.7, 0], 'FontWeight', 'bold');
text(ax_alpha, max(trimmed_time)*0.02, anesthetized_threshold*1.02, 'Anesthetized Threshold', 'Color', [0.7, 0, 0], 'FontWeight', 'bold');
hold(ax_alpha, 'off');

% 3. Window comparison plot at bottom
ax_compare = axes('Parent', tab3, 'Position', [0.05, 0.10, 0.75, 0.25]);
% Initialize with empty plot, will be updated when windows are selected
title(ax_compare, 'Window Comparison: Select two time windows', 'FontWeight', 'bold');
xlabel(ax_compare, 'Frequency (Hz)');
ylabel(ax_compare, 'Power (dB)');
grid(ax_compare, 'on');

% Store the handles
fig.UserData.h_spec3 = h_spec3;
fig.UserData.h_alpha = h_alpha;
fig.UserData.ax_spec3 = ax_spec3;
fig.UserData.ax_alpha = ax_alpha;
fig.UserData.ax_compare = ax_compare;

% Link x-axes for top two plots
linkaxes([ax_spec3, ax_alpha], 'x');

% ===== CONTROLS SECTION (SHARED) =====
control_panel = uipanel('Parent', fig, ...
                        'Title', 'Threshold Controls', ...
                        'Position', [0.86, 0.70, 0.13, 0.25], ... 
                        'FontWeight', 'bold');

colorscale_panel = uipanel('Parent', fig, ...
                        'Title', 'Spectrogram Controls', ...
                        'Position', [0.86, 0.45, 0.13, 0.25], ... 
                        'FontWeight', 'bold');

util_panel = uipanel('Parent', fig, ...
                    'Title', 'Utilities', ...
                    'Position', [0.86, 0.20, 0.13, 0.25], ... 
                    'FontWeight', 'bold');

aperiodic_panel = uipanel('Parent', fig, ...
                    'Title', 'Aperiodic Analysis', ...
                    'Position', [0.86, 0.05, 0.13, 0.15], ... 
                    'FontWeight', 'bold');

% Alpha smoothing slider
uicontrol('Parent', aperiodic_panel, ...
          'Style', 'text', ...
          'String', 'Alpha Smoothing:', ...
          'Position', [10, 105, 100, 15]);

uicontrol('Parent', aperiodic_panel, ...
          'Style', 'slider', ...
          'Min', 1, ...
          'Max', 50, ...
          'Value', 10, ...
          'Position', [10, 90, 150, 15], ...
          'Tag', 'AlphaSmoothingSlider', ...
          'Callback', @(src,event) AlphaSmoothingCallback(src,event));

% Window 1 selection
uicontrol('Parent', aperiodic_panel, ...
          'Style', 'text', ...
          'String', 'Window 1 (min:min):', ...
          'Position', [10, 70, 150, 15]);

uicontrol('Parent', aperiodic_panel, ...
          'Style', 'edit', ...
          'String', [num2str(trimmed_time(1)) ':' num2str(trimmed_time(1)+5)], ...
          'Position', [10, 55, 150, 15], ...
          'Tag', 'Window1Edit');

% Window 2 selection
uicontrol('Parent', aperiodic_panel, ...
          'Style', 'text', ...
          'String', 'Window 2 (min:min):', ...
          'Position', [10, 40, 150, 15]);

uicontrol('Parent', aperiodic_panel, ...
          'Style', 'edit', ...
          'String', [num2str(trimmed_time(end)-5) ':' num2str(trimmed_time(end))], ...
          'Position', [10, 25, 150, 15], ...
          'Tag', 'Window2Edit');

% Compare button
uicontrol('Parent', aperiodic_panel, ...
          'Style', 'pushbutton', ...
          'String', 'Compare Windows', ...
          'Position', [10, 5, 150, 20], ...
          'Callback', @(src,event) CompareWindowsCallback(src,event));

% Threshold type dropdown
uicontrol('Parent', control_panel, ...
          'Style', 'popup', ...
          'String', {'baseline', 'absolute'}, ...
          'Value', strcmp(threshold_type, 'absolute') + 1, ...
          'Position', [10, 175, 120, 25], ...
          'Tag', 'ThresholdTypeDropdown', ...
          'Callback', @(src,event) thresholdTypeCallback(src,event));

% AMI threshold slider and label
uicontrol('Parent', control_panel, ...
          'Style', 'text', ...
          'String', sprintf('AMI Threshold: %.1f', AMI_threshold), ...
          'Position', [10, 150, 120, 20], ...
          'Tag', 'AMIThresholdLabel');

uicontrol('Parent', control_panel, ...
          'Style', 'slider', ...
          'Min', 0, ...
          'Max', 20, ...
          'Value', AMI_threshold, ...
          'Position', [10, 130, 150, 20], ...
          'Tag', 'AMIThresholdSlider', ...
          'Callback', @(src,event) AMIThresholdCallback(src,event));

% SMI threshold slider and label
uicontrol('Parent', control_panel, ...
          'Style', 'text', ...
          'String', sprintf('SMI Threshold: %.1f', SMI_threshold), ...
          'Position', [10, 110, 120, 20], ...
          'Tag', 'SMIThresholdLabel');

uicontrol('Parent', control_panel, ...
          'Style', 'slider', ...
          'Min', 0, ...
          'Max', 20, ...
          'Value', SMI_threshold, ...
          'Position', [10, 90, 150, 20], ...
          'Tag', 'SMIThresholdSlider', ...
          'Callback', @(src,event) SMIThresholdCallback(src,event));

% Add baseline window controls
uicontrol('Parent', control_panel, ...
          'Style', 'text', ...
          'String', 'Baseline Window (min:min):', ...
          'Position', [10, 70, 150, 20]);

% Convert baseline indices to minutes for display
baseline_start_min = trimmed_time(max(1, min(length(trimmed_time), baseline_start_idx)));
baseline_end_min = trimmed_time(max(1, min(length(trimmed_time), baseline_end_idx)));

uicontrol('Parent', control_panel, ...
          'Style', 'edit', ...
          'String', [num2str(baseline_start_min, '%.2f') ':' num2str(baseline_end_min, '%.2f')], ...
          'Position', [10, 50, 150, 20], ...
          'Tag', 'BaselineWindowEdit', ...
          'Callback', @(src,event) BaselineWindowCallback(src,event));

% Apply baseline button
uicontrol('Parent', control_panel, ...
          'Style', 'pushbutton', ...
          'String', 'Apply Baseline Window', ...
          'Position', [10, 30, 150, 20], ...
          'Callback', @(src,event) ApplyBaselineCallback(src,event));


% Min color value slider and label
uicontrol('Parent', colorscale_panel, ...
          'Style', 'text', ...
          'String', sprintf('Min Value (dB): %.1f', cmin), ...
          'Position', [10, 175, 140, 20], ...
          'Tag', 'CMinLabel');

uicontrol('Parent', colorscale_panel, ...
          'Style', 'slider', ...
          'Min', -50, ...
          'Max', 0, ...
          'Value', cmin, ...
          'Position', [10, 155, 150, 20], ...
          'Tag', 'CMinSlider', ...
          'Callback', @(src,event) CMinCallback(src,event));

% Max color value slider and label
uicontrol('Parent', colorscale_panel, ...
          'Style', 'text', ...
          'String', sprintf('Max Value (dB): %.1f', cmax), ...
          'Position', [10, 130, 140, 20], ...
          'Tag', 'CMaxLabel');

uicontrol('Parent', colorscale_panel, ...
          'Style', 'slider', ...
          'Min', 0, ...
          'Max', 50, ...
          'Value', cmax, ...
          'Position', [10, 110, 150, 20], ...
          'Tag', 'CMaxSlider', ...
          'Callback', @(src,event) CMaxCallback(src,event));

% Color map selection with preview
uicontrol('Parent', colorscale_panel, ...
          'Style', 'text', ...
          'String', 'Colormap:', ...
          'Position', [10, 85, 80, 20]);

h_colormap_popup = uicontrol('Parent', colorscale_panel, ...
          'Style', 'popup', ...
          'String', {'jet', 'parula', 'hot', 'cool', 'turbo'}, ...
          'Value', 1, ... % Default to jet
          'Position', [10, 65, 140, 20], ...
          'Callback', @(src,event) ColormapCallback(src,event));

% Add smoothing option for spectrogram
uicontrol('Parent', colorscale_panel, ...
          'Style', 'text', ...
          'String', 'Smoothing:', ...
          'Position', [5, 35, 80, 20]);

uicontrol('Parent', colorscale_panel, ...
          'Style', 'popup', ...
          'String', {'None', 'Light', 'Medium', 'Heavy'}, ...
          'Value', 1, ...
          'Position', [90, 35, 50, 20], ...
          'Callback', @(src,event) SmoothingCallback(src,event));

% Reset zoom button
uicontrol('Parent', util_panel, ...
          'Style', 'pushbutton', ...
          'String', 'Reset Zoom', ...
          'Position', [10, 175, 150, 25], ...
          'Callback', @(src,event) ResetZoomCallback(src,event));

% Export figure button
uicontrol('Parent', util_panel, ...
          'Style', 'pushbutton', ...
          'String', 'Export Figure', ...
          'Position', [10, 145, 150, 25], ...
          'Callback', @(src,event) ExportFigureCallback(src,event));

% Time range selection
uicontrol('Parent', util_panel, ...
          'Style', 'text', ...
          'String', 'Time Range (min:max min):', ...
          'Position', [10, 115, 150, 20]);

uicontrol('Parent', util_panel, ...
          'Style', 'edit', ...
          'String', [num2str(min(trimmed_time)) ':' num2str(max(trimmed_time))], ...
          'Position', [10, 95, 150, 20], ...
          'Callback', @(src,event) TimeRangeCallback(src,event));

% Create section for bandpass controls within the utilities panel
uicontrol('Parent', util_panel, ...
          'Style', 'text', ...
          'String', 'Bandpass Filter:', ...
          'Position', [10, 70, 150, 15], ...
          'FontWeight', 'bold');

% Low and high cutoff
uicontrol('Parent', util_panel, ...
          'Style', 'text', ...
          'String', 'Freq Range (Hz):', ...
          'Position', [5, 50, 90, 15]);

low_cutoff = uicontrol('Parent', util_panel, ...
                       'Style', 'edit', ...
                       'String', '8', ...
                       'Position', [100, 50, 35, 20], ...
                       'Tag', 'LowCutoffEdit');

uicontrol('Parent', util_panel, ...
          'Style', 'text', ...
          'String', 'to', ...
          'Position', [130, 50, 15, 15], ...
          'HorizontalAlignment', 'center');

high_cutoff = uicontrol('Parent', util_panel, ...
                        'Style', 'edit', ...
                        'String', '13', ...
                        'Position', [145, 50, 35, 20], ...
                        'Tag', 'HighCutoffEdit');

% Filter order and type in one row
uicontrol('Parent', util_panel, ...
          'Style', 'text', ...
          'String', 'Order:', ...
          'Position', [10, 25, 40, 15]);

filter_order = uicontrol('Parent', util_panel, ...
                         'Style', 'edit', ...
                         'String', '4', ...
                         'Position', [55, 25, 25, 20], ...
                         'Tag', 'FilterOrderEdit');

uicontrol('Parent', util_panel, ...
          'Style', 'text', ...
          'String', 'Type:', ...
          'Position', [85, 25, 35, 15]);

filter_type = uicontrol('Parent', util_panel, ...
                        'Style', 'popup', ...
                        'String', {'Butter', 'Cheby1', 'Ellip'}, ...
                        'Position', [125, 25, 50, 20], ...
                        'Value', 1, ...
                        'Tag', 'FilterTypeDropdown');

% Apply and Reset buttons on same row
apply_filter = uicontrol('Parent', util_panel, ...
                         'Style', 'pushbutton', ...
                         'String', 'Apply', ...
                         'Position', [10, 5, 80, 20], ...
                         'Callback', @(src,event) ApplyFilterCallback(src,event));

reset_filter = uicontrol('Parent', util_panel, ...
                         'Style', 'pushbutton', ...
                         'String', 'Reset', ...
                         'Position', [100, 5, 80, 20], ...
                         'Callback', @(src,event) ResetFilterCallback(src,event));

% Add status indicator
status_text = uicontrol('Parent', fig, ...
                        'Style', 'text', ...
                        'String', 'Ready', ...
                        'Position', [10, 5, 300, 20]);
fig.UserData.status_text = status_text;

% Display completion message
disp('Enhanced interactive visualization complete. Use the controls to adjust thresholds and spectrogram appearance.');
updateStatus('Visualization ready. Use controls to customize display.', fig);

% CALLBACK FUNCTIONS
    % Callback function for threshold type change
    function thresholdTypeCallback(source, ~)
        fig = ancestor(source, 'figure');
        updateStatus('Processing threshold change...', fig);
        
        val = source.Value;
        threshold_types = {'baseline', 'absolute'};
        selected_type = threshold_types{val};
        
        % Update stored threshold type
        fig.UserData.current_threshold_type = selected_type;
        
        % Clear any existing baseline markers first
        if isfield(fig.UserData, 'baseline_lines') && ~isempty(fig.UserData.baseline_lines)
            for i = 1:length(fig.UserData.baseline_lines)
                if isgraphics(fig.UserData.baseline_lines{i})
                    delete(fig.UserData.baseline_lines{i});
                end
            end
            fig.UserData.baseline_lines = {};
        end
        
        % Recalculate AMI and SMI based on new threshold type
        if strcmp(selected_type, 'baseline')
            [AMI, SMI] = computeAMI_SMI(fig.UserData.eeg_data, fig.UserData.fs, ...
                                      fig.UserData.baseline_start_idx, fig.UserData.baseline_end_idx);
            
            % Get baseline time values for markers
            baseline_start_time = fig.UserData.trimmed_time(max(1, min(length(fig.UserData.trimmed_time), fig.UserData.baseline_start_idx)));
            baseline_end_time = fig.UserData.trimmed_time(max(1, min(length(fig.UserData.trimmed_time), fig.UserData.baseline_end_idx)));
            
            % Add baseline markers to spectrograms
            addBaselineMarkers(fig, baseline_start_time, baseline_end_time);
        else
            [AMI, SMI] = computeAMI_SMI(fig.UserData.eeg_data, fig.UserData.fs, ...
                                      fig.UserData.baseline_start_idx, fig.UserData.baseline_end_idx, ...
                                      'AlphaThresholdType', 'absolute', ...
                                      'AlphaThresholdValue', fig.UserData.AMI_threshold, ...
                                      'SMIThresholdType', 'absolute', ...
                                      'SMIThresholdValue', fig.UserData.SMI_threshold);
        end
        
        % Update plots
        fig.UserData.ami_plot.YData = AMI;
        fig.UserData.smi_plot.YData = SMI;
        
        % Store updated values
        fig.UserData.AMI = AMI;
        fig.UserData.SMI = SMI;
        
        % Enable/disable sliders based on threshold type
        sliders_enabled = strcmp(selected_type, 'absolute');
        ami_slider = findobj(fig, 'Tag', 'AMIThresholdSlider');
        smi_slider = findobj(fig, 'Tag', 'SMIThresholdSlider');
        
        if ~isempty(ami_slider)
            ami_slider.Enable = matlab.lang.OnOffSwitchState(sliders_enabled);
        end
        
        if ~isempty(smi_slider)
            smi_slider.Enable = matlab.lang.OnOffSwitchState(sliders_enabled);
        end
        
        updateStatus('Threshold type updated to: ' + string(selected_type), fig);
    end

    % Callback function for AMI threshold change
    function AMIThresholdCallback(source, ~)
        fig = ancestor(source, 'figure');
        new_threshold = source.Value;
        
        % Update the label
        ami_label = findobj(fig, 'Tag', 'AMIThresholdLabel');
        if ~isempty(ami_label)
            ami_label.String = sprintf('AMI Threshold: %.1f', new_threshold);
        end
        
        % Update stored threshold
        fig.UserData.AMI_threshold = new_threshold;
        
        % Only recalculate if in absolute mode
        if strcmp(fig.UserData.current_threshold_type, 'absolute')
            updateStatus('Recalculating AMI with new threshold...', fig);
            [AMI, ~] = computeAMI_SMI(fig.UserData.eeg_data, fig.UserData.fs, ...
                                    fig.UserData.baseline_start_idx, fig.UserData.baseline_end_idx, ...
                                    'AlphaThresholdType', 'absolute', ...
                                    'AlphaThresholdValue', new_threshold, ...
                                    'SMIThresholdType', 'absolute', ...
                                    'SMIThresholdValue', fig.UserData.SMI_threshold);
            
            % Update AMI plot
            fig.UserData.ami_plot.YData = AMI;
            
            % Store updated values
            fig.UserData.AMI = AMI;
            updateStatus('AMI threshold updated to: ' + string(new_threshold), fig);
        end
    end

    % Callback function for SMI threshold change
    function SMIThresholdCallback(source, ~)
        fig = ancestor(source, 'figure');
        new_threshold = source.Value;
        
        % Update the label
        smi_label = findobj(fig, 'Tag', 'SMIThresholdLabel');
        if ~isempty(smi_label)
            smi_label.String = sprintf('SMI Threshold: %.1f', new_threshold);
        end
        
        % Update stored threshold
        fig.UserData.SMI_threshold = new_threshold;
        
        % Only recalculate if in absolute mode
        if strcmp(fig.UserData.current_threshold_type, 'absolute')
            updateStatus('Recalculating SMI with new threshold...', fig);
            [~, SMI] = computeAMI_SMI(fig.UserData.eeg_data, fig.UserData.fs, ...
                                    fig.UserData.baseline_start_idx, fig.UserData.baseline_end_idx, ...
                                    'AlphaThresholdType', 'absolute', ...
                                    'AlphaThresholdValue', fig.UserData.AMI_threshold, ...
                                    'SMIThresholdType', 'absolute', ...
                                    'SMIThresholdValue', new_threshold);
            
            % Update SMI plot
            fig.UserData.smi_plot.YData = SMI;
            
            % Store updated values
            fig.UserData.SMI = SMI;
            updateStatus('SMI threshold updated to: ' + string(new_threshold), fig);
        end
    end

    
    % Callback for baseline window change
    function BaselineWindowCallback(source, ~)
        % This function validates the baseline window input but doesn't apply it yet
        fig = ancestor(source, 'figure');
        window_str = source.String;
        
        try
            % Parse the time range string (format: min:max)
            parts = strsplit(window_str, ':');
            if length(parts) ~= 2
                error('Invalid format. Use start:end');
            end
            
            start_time = str2double(parts{1});
            end_time = str2double(parts{2});
            
            if isnan(start_time) || isnan(end_time)
                error('Time values must be numbers');
            end
            
            % Validate the range
            all_time = fig.UserData.trimmed_time;
            abs_min = min(all_time);
            abs_max = max(all_time);
            
            if start_time < abs_min
                start_time = abs_min;
            end
            
            if end_time > abs_max
                end_time = abs_max;
            end
            
            if start_time >= end_time
                error('Start time must be less than end time');
            end
            
            % Format and update the text field with corrected values
            source.String = [num2str(start_time, '%.2f') ':' num2str(end_time, '%.2f')];
            
        catch err
            % Reset to current baseline on error
            baseline_start_min = fig.UserData.trimmed_time(max(1, min(length(fig.UserData.trimmed_time), fig.UserData.baseline_start_idx)));
            baseline_end_min = fig.UserData.trimmed_time(max(1, min(length(fig.UserData.trimmed_time), fig.UserData.baseline_end_idx)));
            source.String = [num2str(baseline_start_min, '%.2f') ':' num2str(baseline_end_min, '%.2f')];
            updateStatus(['Error: ' err.message ' - Reset to current baseline'], fig);
        end
    end
    
    % Callback for applying the baseline window
    function ApplyBaselineCallback(~, ~)
        fig = gcf;
        updateStatus('Applying new baseline window...', fig);
        
        % Get the baseline window string
        baseline_edit = findobj(fig, 'Tag', 'BaselineWindowEdit');
        if isempty(baseline_edit)
            updateStatus('Error: Baseline window control not found', fig);
            return;
        end
        
        window_str = baseline_edit.String;
        
        try
            % Parse the time range string
            parts = strsplit(window_str, ':');
            start_time = str2double(parts{1});
            end_time = str2double(parts{2});
            
            % Convert time in minutes to indices
            all_time = fig.UserData.trimmed_time;
            
            % Find closest time points
            [~, start_idx] = min(abs(all_time - start_time));
            [~, end_idx] = min(abs(all_time - end_time));
            
            % Print debug info
            disp(['DEBUG - Selected baseline range: ' num2str(start_time) ' to ' num2str(end_time) ' minutes']);
            disp(['DEBUG - Converted to indices: ' num2str(start_idx) ' to ' num2str(end_idx)]);
            disp(['DEBUG - Corresponding actual times: ' num2str(all_time(start_idx)) ' to ' num2str(all_time(end_idx)) ' minutes']);
            
            % Update stored baseline indices
            fig.UserData.baseline_start_idx = start_idx;
            fig.UserData.baseline_end_idx = end_idx;
            
            % Convert to seconds for AMI/SMI calculation
            baseline_start_sec = start_idx * (length(fig.UserData.eeg_data) / length(all_time) / fig.UserData.fs);
            baseline_end_sec = end_idx * (length(fig.UserData.eeg_data) / length(all_time) / fig.UserData.fs);
            
            % If in baseline threshold mode, recalculate indices
            if strcmp(fig.UserData.current_threshold_type, 'baseline')
                updateStatus('Recalculating indices with new baseline...', fig);
                [AMI, SMI] = computeAMI_SMI(fig.UserData.eeg_data, fig.UserData.fs, ...
                                          baseline_start_sec, baseline_end_sec);
                                      
                % Update plots
                fig.UserData.ami_plot.YData = AMI;
                fig.UserData.smi_plot.YData = SMI;
                
                % Store updated values
                fig.UserData.AMI = AMI;
                fig.UserData.SMI = SMI;
            end
            
            updateStatus(['Baseline window updated to ' window_str ' minutes'], fig);
            
        catch err
            updateStatus(['Error applying baseline: ' err.message], fig);
        end
    end

    % Callback function for minimum colormap value
    function CMinCallback(source, ~)
        fig = ancestor(source, 'figure');
        new_cmin = source.Value;
        
        % Update the label
        cmin_label = findobj(fig, 'Tag', 'CMinLabel');
        if ~isempty(cmin_label)
            cmin_label.String = sprintf('Min Value (dB): %.1f', new_cmin);
        end
        
        % Get current cmax
        cmax = fig.UserData.cmax;
        
        % Ensure cmin is less than cmax
        if new_cmin >= cmax
            new_cmin = cmax - 1;
            source.Value = new_cmin;
            if ~isempty(cmin_label)
                cmin_label.String = sprintf('Min Value (dB): %.1f', new_cmin);
            end
        end
        
        % Update stored value
        fig.UserData.cmin = new_cmin;
        
        % Update the spectrogram colormap limits on both tabs
        set(fig.UserData.ax1, 'clim', [new_cmin, cmax]);
        set(fig.UserData.ax_spec2, 'clim', [new_cmin, cmax]);
        
        updateStatus('Color minimum updated to: ' + string(new_cmin) + ' dB', fig);
    end

    % Callback function for maximum colormap value
    function CMaxCallback(source, ~)
        fig = ancestor(source, 'figure');
        new_cmax = source.Value;
        
        % Update the label
        cmax_label = findobj(fig, 'Tag', 'CMaxLabel');
        if ~isempty(cmax_label)
            cmax_label.String = sprintf('Max Value (dB): %.1f', new_cmax);
        end
        
        % Get current cmin
        cmin = fig.UserData.cmin;
        
        % Ensure cmax is greater than cmin
        if new_cmax <= cmin
            new_cmax = cmin + 1;
            source.Value = new_cmax;
            if ~isempty(cmax_label)
                cmax_label.String = sprintf('Max Value (dB): %.1f', new_cmax);
            end
        end
        
        % Update stored value
        fig.UserData.cmax = new_cmax;
        
        % Update the spectrogram colormap limits on both tabs
        set(fig.UserData.ax1, 'clim', [cmin, new_cmax]);
        set(fig.UserData.ax_spec2, 'clim', [cmin, new_cmax]);
        
        updateStatus('Color maximum updated to: ' + string(new_cmax) + ' dB', fig);
    end

    % Callback function for colormap selection
    function ColormapCallback(source, ~)
        fig = ancestor(source, 'figure');
        val = source.Value;
        colormap_options = {'jet', 'parula', 'hot', 'cool', 'turbo'};
        selected_colormap = colormap_options{val};
        
        % Apply the new colormap to both tabs
        colormap(fig.UserData.ax1, selected_colormap);
        colormap(fig.UserData.ax_spec2, selected_colormap);
        
        updateStatus('Colormap changed to: ' + string(selected_colormap), fig);
    end
    
    % Callback for smoothing
    function SmoothingCallback(source, ~)
        fig = ancestor(source, 'figure');
        smoothing_level = source.Value;
        smoothing_options = {'None', 'Light', 'Medium', 'Heavy'};
        selected_smoothing = smoothing_options{smoothing_level};
        
        updateStatus('Applying ' + string(selected_smoothing) + ' smoothing...', fig);
        
        % Get the original spectrogram data
        spect_data = fig.UserData.spectrogram_data;
        
        % Apply smoothing based on selection
        switch smoothing_level
            case 1 % None
                smoothed_data = spect_data;
            case 2 % Light
                smoothed_data = smooth2d(spect_data, 3);
            case 3 % Medium
                smoothed_data = smooth2d(spect_data, 5);
            case 4 % Heavy
                smoothed_data = smooth2d(spect_data, 9);
        end
        
        % Update both spectrograms
        fig.UserData.h_spec.CData = pow2db(smoothed_data);
        fig.UserData.h_spec2.CData = pow2db(smoothed_data);
        
        updateStatus(selected_smoothing + ' smoothing applied', fig);
    end
    
    % Callback for reset zoom button
    function ResetZoomCallback(~, ~)
        fig = gcf;
        updateStatus('Resetting zoom...', fig);
        
        % Reset zoom for all axes in both tabs
        axes_handles = [fig.UserData.ax1, fig.UserData.ax2, fig.UserData.ax3, ...
                       fig.UserData.ax_eeg, fig.UserData.ax_spec2];
        
        for i = 1:length(axes_handles)
            if isgraphics(axes_handles(i))
                axes(axes_handles(i));
                axis auto;
            end
        end
        
        % Reset time range to full range
        min_time = min(fig.UserData.trimmed_time);
        max_time = max(fig.UserData.trimmed_time);
        
        for i = 1:length(axes_handles)
            if isgraphics(axes_handles(i))
                xlim(axes_handles(i), [min_time, max_time]);
            end
        end
        
        % Update time range display
        time_range_edit = findobj(fig, 'Style', 'edit');
        if ~isempty(time_range_edit)
            time_range_edit.String = [num2str(min_time) ':' num2str(max_time)];
        end
        
        updateStatus('Zoom reset complete', fig);
    end
    
    % Callback for export figure button
    function ExportFigureCallback(~, ~)
        fig = gcf;
        updateStatus('Preparing figure export...', fig);
        
        % Create a dialog to select export options
        export_dialog = figure('Name', 'Export Options', 'Position', [300, 300, 400, 200], ...
                         'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', ...
                         'WindowStyle', 'modal');
        
        % File format selection
        uicontrol('Parent', export_dialog, 'Style', 'text', 'String', 'File Format:', ...
                 'Position', [20, 160, 120, 20], 'HorizontalAlignment', 'left');
        format_popup = uicontrol('Parent', export_dialog, 'Style', 'popup', ...
                           'String', {'PNG (.png)', 'TIFF (.tif)', 'JPEG (.jpg)', 'PDF (.pdf)', 'EPS (.eps)', 'FIG (.fig)'}, ...
                           'Position', [150, 160, 150, 25], 'Value', 1);
        
        % Resolution selection
        uicontrol('Parent', export_dialog, 'Style', 'text', 'String', 'Resolution (DPI):', ...
                 'Position', [20, 120, 120, 20], 'HorizontalAlignment', 'left');
        resolution_popup = uicontrol('Parent', export_dialog, 'Style', 'popup', ...
                               'String', {'72 (Screen)', '150 (Medium)', '300 (Print)', '600 (High)'}, ...
                               'Position', [150, 120, 150, 25], 'Value', 3);
        
        % File name input
        uicontrol('Parent', export_dialog, 'Style', 'text', 'String', 'File Name:', ...
                 'Position', [20, 80, 120, 20], 'HorizontalAlignment', 'left');
        filename_edit = uicontrol('Parent', export_dialog, 'Style', 'edit', ...
                            'String', 'EEG_Analysis_Export', 'Position', [150, 80, 200, 25]);
        
        % Export and Cancel buttons
        uicontrol('Parent', export_dialog, 'Style', 'pushbutton', 'String', 'Export', ...
                 'Position', [100, 30, 80, 30], 'Callback', @(src,event) export_figure_callback(src, event));
        uicontrol('Parent', export_dialog, 'Style', 'pushbutton', 'String', 'Cancel', ...
                 'Position', [220, 30, 80, 30], 'Callback', @(src,event) delete(export_dialog));
        
        % Nested export callback
        function export_figure_callback(~, ~)
            % Get export options
            format_val = format_popup.Value;
            format_options = {'.png', '.tif', '.jpg', '.pdf', '.eps', '.fig'};
            selected_format = format_options{format_val};
            
            resolution_val = resolution_popup.Value;
            resolution_options = [72, 150, 300, 600];
            selected_resolution = resolution_options(resolution_val);
            
            filename = filename_edit.String;
            if isempty(filename)
                filename = 'EEG_Analysis_Export';
            end
            
            full_filename = [filename, selected_format];
            
            % Close the dialog
            delete(export_dialog);
            
            % Export the figure
            updateStatus(['Exporting to ' full_filename '...'], fig);
            
            try
                if strcmp(selected_format, '.fig')
                    % For .fig format, save the entire figure
                    savefig(fig, full_filename);
                else
                    % For image formats, use export_fig (if available) or print
                    if exist('export_fig', 'file') == 2
                        export_fig(fig, full_filename, ['-r' num2str(selected_resolution)]);
                    else
                        print(fig, full_filename, ['-d' selected_format(2:end)], ['-r' num2str(selected_resolution)]);
                    end
                end
                updateStatus(['Figure exported successfully to ' full_filename], fig);
            catch err
                updateStatus(['Error exporting figure: ' err.message], fig);
            end
        end
    end
    
    % Callback for time range input
    function TimeRangeCallback(source, ~)
        fig = gcf;
        time_range_str = source.String;
        
        try
            % Parse the time range string (format: min:max)
            parts = strsplit(time_range_str, ':');
            if length(parts) ~= 2
                error('Invalid time range format. Use min:max');
            end
            
            min_time = str2double(parts{1});
            max_time = str2double(parts{2});
            
            if isnan(min_time) || isnan(max_time)
                error('Time values must be numbers');
            end
            
            % Validate the range
            all_time = fig.UserData.trimmed_time;
            abs_min = min(all_time);
            abs_max = max(all_time);
            
            if min_time < abs_min
                min_time = abs_min;
            end
            
            if max_time > abs_max
                max_time = abs_max;
            end
            
            if min_time >= max_time
                error('Minimum time must be less than maximum time');
            end
            
            % Update the display range for all plots
            updateStatus('Updating time range...', fig);
            
            % Get all axes handles
            axes_handles = [fig.UserData.ax1, fig.UserData.ax2, fig.UserData.ax3, ...
                           fig.UserData.ax_eeg, fig.UserData.ax_spec2];
            
            for i = 1:length(axes_handles)
                if isgraphics(axes_handles(i))
                    xlim(axes_handles(i), [min_time, max_time]);
                end
            end
            
            % Update the text field with corrected values
            source.String = [num2str(min_time) ':' num2str(max_time)];
            updateStatus(['Time range set to ' source.String], fig);
            
        catch err
            % Reset to full range on error
            all_time = fig.UserData.trimmed_time;
            source.String = [num2str(min(all_time)) ':' num2str(max(all_time))];
            updateStatus(['Error: ' err.message ' - Reset to full range'], fig);
        end
    end
    
    % Callback for adding event markers
    function AddMarkerCallback(~, ~)
        fig = gcf;
        updateStatus('Click on a plot to add an event marker...', fig);
        
        % Store the original button down function
        fig.UserData.original_buttondown = fig.WindowButtonDownFcn;
        
        % Set the button down function to handle marker placement
        fig.WindowButtonDownFcn = @PlaceMarker;
        
        % Nested function to handle marker placement
        function PlaceMarker(~, ~)
            % Get the clicked axes and position
            clicked_obj = gca;
            
            % Get the axes position
            point = clicked_obj.CurrentPoint;
            x_pos = point(1, 1);
            
            % Add marker to event list
            if ~isfield(fig.UserData, 'event_markers') || isempty(fig.UserData.event_markers)
                fig.UserData.event_markers = x_pos;
                marker_num = 1;
            else
                fig.UserData.event_markers(end+1) = x_pos;
                marker_num = length(fig.UserData.event_markers);
            end
            
            % Draw marker lines on all plots
            axes_handles = [fig.UserData.ax1, fig.UserData.ax2, fig.UserData.ax3, ...
                           fig.UserData.ax_eeg, fig.UserData.ax_spec2];
            
            for i = 1:length(axes_handles)
                if isgraphics(axes_handles(i))
                    % Create the line
                    hold(axes_handles(i), 'on');
                    marker_line = line(axes_handles(i), [x_pos, x_pos], get(axes_handles(i), 'YLim'), ...
                                      'Color', [1, 0, 0, 0.7], 'LineWidth', 1.5, 'LineStyle', '-', ...
                                      'DisplayName', ['Marker ' num2str(marker_num)]);
                    
                    % Add to the list of marker lines
                    if ~isfield(fig.UserData, 'marker_lines') || isempty(fig.UserData.marker_lines)
                        fig.UserData.marker_lines = {marker_line};
                    else
                        fig.UserData.marker_lines{end+1} = marker_line;
                    end
                    
                    % Add marker label (only for the top plot in each tab)
                    if i == 1 || i == 4
                        y_lim = get(axes_handles(i), 'YLim');
                        y_pos = y_lim(2) - 0.05 * (y_lim(2) - y_lim(1));
                        
                        marker_text = text(axes_handles(i), x_pos, y_pos, ['M' num2str(marker_num)], ...
                                         'Color', 'red', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                        
                        fig.UserData.marker_lines{end+1} = marker_text;
                    end
                    
                    hold(axes_handles(i), 'off');
                end
            end
            
            % Restore original button down function
            fig.WindowButtonDownFcn = fig.UserData.original_buttondown;
            
            updateStatus(['Event marker ' num2str(marker_num) ' added at time ' num2str(x_pos) ' min'], fig);
        end
    end
    
    % Callback for clearing all markers
    function ClearMarkersCallback(~, ~)
        fig = gcf;
        
        % Delete all marker lines
        if isfield(fig.UserData, 'marker_lines') && ~isempty(fig.UserData.marker_lines)
            for i = 1:length(fig.UserData.marker_lines)
                if isgraphics(fig.UserData.marker_lines{i})
                    delete(fig.UserData.marker_lines{i});
                end
            end
            fig.UserData.marker_lines = {};
        end
        
        % Clear marker positions
        fig.UserData.event_markers = [];
        
        updateStatus('All event markers cleared', fig);
    end
    
    % Callback function for applying bandpass filter
function ApplyFilterCallback(~, ~)
    fig = gcf;
    updateStatus('Applying bandpass filter...', fig);
    
    try
        % Get filter parameters from UI controls
        low_freq = str2double(get(findobj(fig, 'Tag', 'LowCutoffEdit'), 'String'));
        high_freq = str2double(get(findobj(fig, 'Tag', 'HighCutoffEdit'), 'String'));
        order = str2double(get(findobj(fig, 'Tag', 'FilterOrderEdit'), 'String'));
        
        % Get filter type
        filter_type_control = findobj(fig, 'Tag', 'FilterTypeDropdown');
        filter_types = {'butter', 'cheby1', 'ellip'};
        selected_filter = filter_types{filter_type_control.Value};
        
        % Validate parameters
        if isnan(low_freq) || isnan(high_freq) || isnan(order)
            error('Filter parameters must be numeric');
        end
        
        if low_freq >= high_freq
            error('Low cutoff must be less than high cutoff');
        end
        
        if low_freq <= 0
            error('Low cutoff must be greater than 0');
        end
        
        % Get sampling frequency and data
        fs = fig.UserData.fs;
        eeg_data = fig.UserData.original_eeg_data;
        
        % Normalize frequencies to Nyquist
        Wn = [low_freq, high_freq] / (fs/2);
        
        % Design the filter
        switch selected_filter
            case 'butter'
                [b, a] = butter(order, Wn, 'bandpass');
            case 'cheby1'
                [b, a] = cheby1(order, 0.5, Wn, 'bandpass');
            case 'ellip'
                [b, a] = ellip(order, 0.5, 40, Wn, 'bandpass');
        end
        
        % Apply the filter
        filtered_data = filtfilt(b, a, double(eeg_data));
        
        % Update the EEG plot with filtered data
        fig.UserData.raw_eeg_plot.YData = filtered_data;
        
        % Store filtered data
        fig.UserData.filtered_data = filtered_data;
        
        % Calculate power in the selected frequency band over time
        % Use a sliding window approach
        window_size = round(fs * 2);  % 2-second window
        overlap = round(window_size * 0.8);  % 80% overlap
        step_size = window_size - overlap;
        
        % Initialize power array
        num_steps = floor((length(filtered_data) - window_size) / step_size) + 1;
        power_values = zeros(1, num_steps);
        time_points = zeros(1, num_steps);
        
        % Calculate power in each window
        for i = 1:num_steps
            start_idx = (i-1) * step_size + 1;
            end_idx = start_idx + window_size - 1;
            
            if end_idx > length(filtered_data)
                end_idx = length(filtered_data);
            end
            
            % Get data segment
            segment = filtered_data(start_idx:end_idx);
            
            % Calculate power (mean squared amplitude)
            power_values(i) = mean(segment.^2);
            
            % Calculate corresponding time point (center of window)
            window_center_idx = start_idx + floor(window_size/2);
            if window_center_idx <= length(fig.UserData.trimmed_time)
                time_points(i) = fig.UserData.trimmed_time(window_center_idx);
            else
                time_points(i) = fig.UserData.trimmed_time(end);
            end
        end
        
        % Apply multi-pass smoothing of the power trajectory
        % First use a wide moving average to get the general trend
        power_values = smoothdata(power_values, 'movmean', 30);
        % Then apply Gaussian smoothing for a more natural appearance
        power_values = smoothdata(power_values, 'gaussian', 15);
        % Finally a little additional smoothing to remove any remaining artifacts
        power_values = smoothdata(power_values, 'lowess', 10);

        % Update power trajectory plot
        fig.UserData.power_plot.XData = time_points;
        fig.UserData.power_plot.YData = power_values;

        % Store band power data
        fig.UserData.bandpower_data = power_values;

        % Update title of the power plot
        title(fig.UserData.ax_power, sprintf('Power in %.1f-%.1f Hz Band (Smoothed)', low_freq, high_freq), 'FontWeight', 'bold');

        % Update y-axis limits for better visualization
        ylim(fig.UserData.ax_power, [0, max(power_values)*1.1]);
        
        % Update spectrogram with new data
        % Recalculate multitaper spectrogram on filtered data
        win_length = 2; % From original code
        nw = win_length * fs;
        N = floor(length(filtered_data) / nw);
        
        if N > 0
            % Reshape filtered data
            yy_filtered = reshape(filtered_data(1:nw*N), nw, N);
            
            % Time-bandwidth product and number of tapers
            TW = 1;
            K = 3;
            
            % Get the tapers
            [tapers, ~] = dpss(nw, TW, K);
            
            % Initialize arrays
            mtY = zeros(nw, K, N);
            for i = 1:N
                for k = 1:K
                    mtY(:, k, i) = tapers(:, k) .* yy_filtered(:, i);
                end
            end
            
            % Compute FFT
            mtFreqY = fft(mtY, nw, 1);
            
            % Average over tapers
            spectrogram_filtered = zeros(nw, N);
            for i = 1:N
                spectrogram_filtered(:, i) = mean(abs(mtFreqY(:, :, i)).^2, 2);
            end
            
            % Update spectrogram plot
            fmax = 30;  % From original code
            fig.UserData.h_spec2.CData = pow2db(spectrogram_filtered(1:fmax*win_length, :));
        end
        
        updateStatus(['Bandpass filter applied: ' num2str(low_freq) '-' num2str(high_freq) ' Hz'], fig);
        
    catch err
        updateStatus(['Error applying filter: ' err.message], fig);
    end
end

% Callback function for resetting to original data
function ResetFilterCallback(~, ~)
    fig = gcf;
    updateStatus('Resetting to original data...', fig);
    
    % Reset EEG plot to original data
    fig.UserData.raw_eeg_plot.YData = fig.UserData.original_eeg_data;
    
    % Reset spectrogram to original
    fig.UserData.h_spec2.CData = pow2db(fig.UserData.original_spectrogram);
    
    % Clear power plot
    fig.UserData.power_plot.YData = zeros(size(fig.UserData.trimmed_time));
    
    % Reset title
    title(fig.UserData.ax_power, 'Bandpass Filtered Power', 'FontWeight', 'bold');
    
    updateStatus('Reset to original data complete', fig);
end

% Callback for alpha smoothing
function AlphaSmoothingCallback(source, ~)
    fig = ancestor(source, 'figure');
    window_size = round(source.Value);
    
    % Get original spectral exponent
    spectral_exponent = fig.UserData.aperiodic_data.spectral_exponent;
    
    % Apply smoothing
    smoothed_alpha = smoothdata(spectral_exponent, 'movmean', window_size);
    
    % Update the plot
    fig.UserData.h_alpha.YData = smoothed_alpha;
    
    updateStatus(['Alpha smoothing window set to ' num2str(window_size)], fig);
end

% Callback for window comparison
function CompareWindowsCallback(~, ~)
    fig = gcf;
    updateStatus('Comparing time windows...', fig);
    
    try
        % Get window ranges
        window1_edit = findobj(fig, 'Tag', 'Window1Edit');
        window2_edit = findobj(fig, 'Tag', 'Window2Edit');
        
        window1_str = window1_edit.String;
        window2_str = window2_edit.String;
        
        % Parse time ranges
        window1_range = str2num(window1_str);
        window2_range = str2num(window2_str);

        window1_range = sscanf(window1_str,'%d:%d');
        window2_range = sscanf(window2_str,'%d:%d');
        
        if length(window1_range) ~= 2 || length(window2_range) ~= 2
            error('Time ranges must be in format start:end');
        end
        
        % Get aperiodic data
        aperiodic_data = fig.UserData.aperiodic_data;
        exponent_time = aperiodic_data.exponent_time;  
        spectral_exponent = aperiodic_data.spectral_exponent;
        aperiodic_spectra = aperiodic_data.aperiodic_spectra;
        freqs = aperiodic_data.freqs;
        
        % Find indices for each window
        window1_idx = exponent_time >= window1_range(1) & exponent_time <= window1_range(2);
        window2_idx = exponent_time >= window2_range(1) & exponent_time <= window2_range(2);
        
        window1_idxx = trimmed_time >= window1_range(1) & trimmed_time <= window1_range(2);
        window2_idxx = trimmed_time >= window2_range(1) & trimmed_time <= window2_range(2);

        if sum(window1_idx) < 1 || sum(window2_idx) < 1
            error('Selected time windows contain no data points');
        end
        
        % Calculate average aperiodic spectra for each window
        avg_spectra1 = mean(aperiodic_spectra(:, window1_idx), 2);
        avg_spectra2 = mean(aperiodic_spectra(:, window2_idx), 2);
        
        % Calculate average spectral exponents
        avg_alpha1 = mean(spectral_exponent(window1_idx));
        avg_alpha2 = mean(spectral_exponent(window2_idx));
        
        % Plot on the comparison axes
        ax_compare = fig.UserData.ax_compare;
        cla(ax_compare);
        hold(ax_compare, 'on');
        
        % Plot average spectra
        plot(ax_compare, freqs, pow2db(avg_spectra1), 'b-', 'LineWidth', 2);
        plot(ax_compare, freqs, pow2db(avg_spectra2), 'r-', 'LineWidth', 2);

        [p1,f1] = pmtm(eeg_data(window1_idxx),10,700,fs);
        [p2,f2] = pmtm(eeg_data(window2_idxx),10,700,fs);

        % p1 = abs(fft(eeg_data(window1_idx))); N = numel(eeg_data(window1_idx)); f1 = fs*(0:(N-1/fs))/N;
        % [p2,f2] = pspectrum(eeg_data(window2_idx),fs,'FrequencyResolution',2);

        plot(ax_compare,f1,pow2db(p1),'-r');
        plot(ax_compare,f2,pow2db(p2),'-b');

        xlim([0,40])
        % plot(ax_compare,freqs, pow2db(pspectrum(eeg_data(window1_idx),s, 'FrequencyLimits',[freqs(1) freqs(end)])),'-r');
        % plot(ax_compare,freqs, (pwelch(eeg_data(window2_idx),[],[],freqs,fs)),'-b');

        % Add legend with alpha values
        legend(ax_compare, ...
               {['Window 1: \alpha = ' num2str(avg_alpha1, '%.2f')], ...
                ['Window 2: \alpha = ' num2str(avg_alpha2, '%.2f')]}, ...
               'Location', 'northeast');
        
        % Update title
        title(ax_compare, sprintf('Aperiodic Spectra Comparison (\\Delta\\alpha = %.2f)', avg_alpha2 - avg_alpha1), ...
              'FontWeight', 'bold');
        
        % Format axes
        xlabel(ax_compare, 'Frequency (Hz)');
        ylabel(ax_compare, 'Power (dB)');
        grid(ax_compare, 'on');
        hold(ax_compare, 'off');
        
        % Add markers on the alpha plot
        ax_alpha = fig.UserData.ax_alpha;
        
        % Clear any existing window markers
        lines_to_keep = findobj(ax_alpha, 'Type', 'line', 'LineStyle', '--', 'Color', [0, 0.7, 0]);
        lines_to_keep = [lines_to_keep; findobj(ax_alpha, 'Type', 'line', 'LineStyle', '--', 'Color', [0.7, 0, 0])];
        delete(setdiff(findobj(ax_alpha, 'Type', 'line'), [lines_to_keep; fig.UserData.h_alpha]));
        
        % Add new window markers
        hold(ax_alpha, 'on');
        y_limits = get(ax_alpha, 'YLim');
        
        % Window 1 boundaries (blue)
        plot(ax_alpha, [window1_range(1), window1_range(1)], y_limits, 'b-', 'LineWidth', 1.5);
        plot(ax_alpha, [window1_range(2), window1_range(2)], y_limits, 'b-', 'LineWidth', 1.5);
        
        % Window 2 boundaries (red)
        plot(ax_alpha, [window2_range(1), window2_range(1)], y_limits, 'r-', 'LineWidth', 1.5);
        plot(ax_alpha, [window2_range(2), window2_range(2)], y_limits, 'r-', 'LineWidth', 1.5);
        
        hold(ax_alpha, 'off');
        
        updateStatus('Window comparison complete', fig);
        
    catch err
        updateStatus(['Error comparing windows: ' err.message], fig);
    end
end

% Function to update status message
    function updateStatus(message, fig)
        if isfield(fig.UserData, 'status_text') && isgraphics(fig.UserData.status_text)
            fig.UserData.status_text.String = message;
            drawnow;
        end
    end
    
    % 2D smoothing function for spectrograms
    function smoothed = smooth2d(data, window_size)
        % Create a 2D averaging filter
        filter = ones(window_size) / window_size^2;
        
        % Apply the filter using conv2
        smoothed = conv2(data, filter, 'same');
    end
end