% Load data from the CSV file with original headers preserved
data = readtable('prop_conc_data_new.csv', 'VariableNamingRule', 'preserve');

% Extract relevant columns using original headers
time = data.("time(sec)");
blood_concentration = data.("Conc blood(C1)")/1000;
effect_site_concentration = data.("Conc effect")/1000;
% predicted_BIS = data.("Predicted BIS");

% Create a figure
figure;

% Create the left y-axis plot (Blood and Effect-site Concentrations)
% yyaxis left;
h1 = plot((time/60), blood_concentration, 'r-', 'LineWidth', 2); % Blood Concentration (Red)
hold on;
h2 = plot((time/60), effect_site_concentration, 'b--', 'LineWidth', 2); % Effect-site Concentration (Blue)
ylabel('Concentration (mcg/ml)', 'FontSize', 15);
xlabel('Time (minutes)', 'FontSize', 15);
grid on;

% Create the right y-axis plot (Predicted BIS)
% yyaxis right;
% h3 = plot((time/60), predicted_BIS, 'g-.', 'LineWidth', 2); % Predicted BIS (Green)
% ylabel('Predicted BIS', 'FontSize', 12);

% Add a title
title('Blood concentration and Effect-site concentration', 'FontSize', 14);

% Explicitly define the legend entries
% legend([h1, h2, h3], 'Blood Concentration (Cp)', 'Effect-site Concentration (Ce)', 'Predicted BIS', ...
%     'Location', 'NorthWest');
legend([h1, h2], 'Blood Concentration (Cp)', 'Effect-site Concentration (Ce)', ...
    'Location', 'NorthWest');

% Final touches
set(gca, 'FontSize', 15);
