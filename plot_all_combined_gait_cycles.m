% Process all AMC files and plot all gait cycles together
% This script finds all .amc files, extracts gait cycles, and plots them all overlaid

close all; clear all; clc
%% Get all AMC files in AMC folder
fprintf('Searching for AMC files in AMC folder...\n');
amc_files = dir('AMC/*.amc');
% amc_path = addpath('./AMC/');
if isempty(amc_files)
    error('No AMC files found in current directory!');
end

fprintf('Found %d AMC files:\n', length(amc_files));
for i = 1:length(amc_files)
    fprintf('  %d. %s\n', i, amc_files(i).name);
end


%% Initialize storage for all cycles from all files
all_right_hip_cycles = [];
all_left_hip_cycles = [];
all_right_knee_cycles = [];
all_left_knee_cycles = [];

% Store file index for each cycle to track which file it came from
right_hip_file_indices = [];
left_hip_file_indices = [];
right_knee_file_indices = [];
left_knee_file_indices = [];

% Color scheme for different files
file_colors = lines(length(amc_files));

%% Process each AMC file and collect all cycles
fprintf('\n=== PROCESSING ALL FILES ===\n');
total_cycles_found = 0;

for file_idx = 1:length(amc_files)
    filename = amc_files(file_idx).name;
    fprintf('Processing file %d/%d: %s\n', file_idx, length(amc_files), filename);

    try
        % Extract gait cycles using robust method
        gait_cycles_data = extract_gait_cycles_knee_minima_robust(filename);

        if isempty(gait_cycles_data)
            fprintf('  WARNING: No gait cycles detected in %s\n', filename);
            continue;
        end

        fprintf('  Found %d gait cycles\n', length(gait_cycles_data));
        total_cycles_found = total_cycles_found + length(gait_cycles_data);

        % Separate cycles by leg
        right_cycles = find(strcmp({gait_cycles_data.leg}, 'right'));
        left_cycles = find(strcmp({gait_cycles_data.leg}, 'left'));

        % Interpolation parameters
        interp_length = 200; % 0 to 100% in 1% steps
        time_standard = linspace(0, 1, interp_length);

        % Process right leg cycles
        for i = 1:length(right_cycles)
            cycle_idx = right_cycles(i);
            time_norm = gait_cycles_data(cycle_idx).time_normalized;

            % Interpolate to standard length
            right_hip_interp = interp1(time_norm, gait_cycles_data(cycle_idx).right_hip_flex, time_standard, 'linear');
            right_knee_interp = interp1(time_norm, gait_cycles_data(cycle_idx).right_knee_flex, time_standard, 'linear');

            % Store interpolated data
            all_right_hip_cycles = [all_right_hip_cycles; right_hip_interp];
            all_right_knee_cycles = [all_right_knee_cycles; right_knee_interp];

            % Store file index for this cycle
            right_hip_file_indices = [right_hip_file_indices; file_idx];
            right_knee_file_indices = [right_knee_file_indices; file_idx];
        end

        % Process left leg cycles
        for i = 1:length(left_cycles)
            cycle_idx = left_cycles(i);
            time_norm = gait_cycles_data(cycle_idx).time_normalized;

            % Interpolate to standard length
            left_hip_interp = interp1(time_norm, gait_cycles_data(cycle_idx).left_hip_flex, time_standard, 'linear');
            left_knee_interp = interp1(time_norm, gait_cycles_data(cycle_idx).left_knee_flex, time_standard, 'linear');

            % Store interpolated data
            all_left_hip_cycles = [all_left_hip_cycles; left_hip_interp];
            all_left_knee_cycles = [all_left_knee_cycles; left_knee_interp];

            % Store file index for this cycle
            left_hip_file_indices = [left_hip_file_indices; file_idx];
            left_knee_file_indices = [left_knee_file_indices; file_idx];
        end

    catch ME
        fprintf('  ERROR processing %s: %s\n', filename, ME.message);
    end
end

fprintf('\nTotal cycles collected: %d\n', total_cycles_found);
fprintf('Right hip cycles: %d\n', size(all_right_hip_cycles, 1));
fprintf('Left hip cycles: %d\n', size(all_left_hip_cycles, 1));
fprintf('Right knee cycles: %d\n', size(all_right_knee_cycles, 1));
fprintf('Left knee cycles: %d\n', size(all_left_knee_cycles, 1));

if total_cycles_found == 0
    error('No gait cycles found in any files!');
end

%% Create comprehensive plot with all cycles from all files
figure('Name', 'All Gait Cycles from All AMC Files', 'Position', [0, 50, 1550, 800]);

time_standard = linspace(0, 1, 200);

%% Plot 1: Right Hip Flexion - All cycles from all files
subplot(2,3,1);
hold on;

for i = 1:size(all_right_hip_cycles, 1)
    file_idx = right_hip_file_indices(i);
    plot(time_standard, all_right_hip_cycles(i, :), ...
         'Color', [file_colors(file_idx, :), 0.6], 'LineWidth', 1);
end
xlabel('Normalized Time (0-1)');
ylabel('Right Hip Flexion (degrees)');
title(sprintf('Right Hip Flexion - All %d Cycles', size(all_right_hip_cycles, 1)));
grid on;

% Add legend for files
legend_entries = {};
for file_idx = 1:length(amc_files)
    [~, name_only, ~] = fileparts(amc_files(file_idx).name);
    legend_entries{end+1} = name_only;
    plot(NaN, NaN, 'Color', file_colors(file_idx, :), 'LineWidth', 2);
end
legend(legend_entries, 'Location', 'best', 'FontSize', 8);

%% Plot 2: Left Hip Flexion - All cycles from all files
subplot(2,3,2);
hold on;

for i = 1:size(all_left_hip_cycles, 1)
    file_idx = left_hip_file_indices(i);
    plot(time_standard, all_left_hip_cycles(i, :), ...
         'Color', [file_colors(file_idx, :), 0.6], 'LineWidth', 1);
end
xlabel('Normalized Time (0-1)');
ylabel('Left Hip Flexion (degrees)');
title(sprintf('Left Hip Flexion - All %d Cycles', size(all_left_hip_cycles, 1)));
grid on;

% Add legend for files
for file_idx = 1:length(amc_files)
    plot(NaN, NaN, 'Color', file_colors(file_idx, :), 'LineWidth', 2);
end
legend(legend_entries, 'Location', 'best', 'FontSize', 8);

%% Plot 3: Hip Flexion - Mean + STD for all cycles combined
subplot(2,3,3);
hold on;

if ~isempty(all_right_hip_cycles)
    right_hip_mean = mean(all_right_hip_cycles, 1);
    right_hip_std = std(all_right_hip_cycles, 0, 1);

    fill([time_standard, fliplr(time_standard)], ...
         [right_hip_mean + right_hip_std, fliplr(right_hip_mean - right_hip_std)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Right ± STD');
    plot(time_standard, right_hip_mean, 'r-', 'LineWidth', 3, 'DisplayName', 'Right Mean');
end

if ~isempty(all_left_hip_cycles)
    left_hip_mean = mean(all_left_hip_cycles, 1);
    left_hip_std = std(all_left_hip_cycles, 0, 1);

    fill([time_standard, fliplr(time_standard)], ...
         [left_hip_mean + left_hip_std, fliplr(left_hip_mean - left_hip_std)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Left ± STD');
    plot(time_standard, left_hip_mean, 'b-', 'LineWidth', 3, 'DisplayName', 'Left Mean');
end

xlabel('Normalized Time (0-1)');
ylabel('Hip Flexion (degrees)');
title('Hip Flexion - Combined Mean ± STD');
legend('Location', 'best');
grid on;

%% Plot 4: Right Knee Flexion - All cycles from all files
subplot(2,3,4);
hold on;

for i = 1:size(all_right_knee_cycles, 1)
    file_idx = right_knee_file_indices(i);
    plot(time_standard, all_right_knee_cycles(i, :), ...
         'Color', [file_colors(file_idx, :), 0.6], 'LineWidth', 1);
end
xlabel('Normalized Time (0-1)');
ylabel('Right Knee Flexion (degrees)');
title(sprintf('Right Knee Flexion - All %d Cycles', size(all_right_knee_cycles, 1)));
grid on;

% Add legend for files
for file_idx = 1:length(amc_files)
    plot(NaN, NaN, 'Color', file_colors(file_idx, :), 'LineWidth', 2);
end
legend(legend_entries, 'Location', 'best', 'FontSize', 8);

%% Plot 5: Left Knee Flexion - All cycles from all files
subplot(2,3,5);
hold on;

for i = 1:size(all_left_knee_cycles, 1)
    file_idx = left_knee_file_indices(i);
    plot(time_standard, all_left_knee_cycles(i, :), ...
         'Color', [file_colors(file_idx, :), 0.6], 'LineWidth', 1);
end
xlabel('Normalized Time (0-1)');
ylabel('Left Knee Flexion (degrees)');
title(sprintf('Left Knee Flexion - All %d Cycles', size(all_left_knee_cycles, 1)));
grid on;

% Add legend for files
for file_idx = 1:length(amc_files)
    plot(NaN, NaN, 'Color', file_colors(file_idx, :), 'LineWidth', 2);
end
legend(legend_entries, 'Location', 'best', 'FontSize', 8);

%% Plot 6: Knee Flexion - Mean + STD for all cycles combined
subplot(2,3,6);
hold on;

if ~isempty(all_right_knee_cycles)
    right_knee_mean = mean(all_right_knee_cycles, 1);
    right_knee_std = std(all_right_knee_cycles, 0, 1);

    fill([time_standard, fliplr(time_standard)], ...
         [right_knee_mean + right_knee_std, fliplr(right_knee_mean - right_knee_std)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Right ± STD');
    plot(time_standard, right_knee_mean, 'r-', 'LineWidth', 3, 'DisplayName', 'Right Mean');
end

if ~isempty(all_left_knee_cycles)
    left_knee_mean = mean(all_left_knee_cycles, 1);
    left_knee_std = std(all_left_knee_cycles, 0, 1);

    fill([time_standard, fliplr(time_standard)], ...
         [left_knee_mean + left_knee_std, fliplr(left_knee_mean - left_knee_std)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Left ± STD');
    plot(time_standard, left_knee_mean, 'b-', 'LineWidth', 3, 'DisplayName', 'Left Mean');
end

xlabel('Normalized Time (0-1)');
ylabel('Knee Flexion (degrees)');
title('Knee Flexion - Combined Mean ± STD');
legend('Location', 'best');
grid on;

%% Add overall title
sgtitle(sprintf('Combined Gait Analysis: %d Files, %d Total Cycles', length(amc_files), total_cycles_found));

%% Save the combined plot
combined_fig_filename = sprintf('combined_gait_analysis_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
saveas(gcf, combined_fig_filename);
fprintf('\nCombined plot saved as: %s\n', combined_fig_filename);

%% Save all data to MAT file
combined_data = struct();
combined_data.amc_files = {amc_files.name};
combined_data.all_right_hip_cycles = all_right_hip_cycles;
combined_data.all_left_hip_cycles = all_left_hip_cycles;
combined_data.all_right_knee_cycles = all_right_knee_cycles;
combined_data.all_left_knee_cycles = all_left_knee_cycles;
combined_data.right_hip_file_indices = right_hip_file_indices;
combined_data.left_hip_file_indices = left_hip_file_indices;
combined_data.right_knee_file_indices = right_knee_file_indices;
combined_data.left_knee_file_indices = left_knee_file_indices;
combined_data.time_standard = time_standard;
combined_data.processing_time = datetime('now');

% Calculate combined statistics
if ~isempty(all_right_hip_cycles)
    combined_data.right_hip_mean = mean(all_right_hip_cycles, 1);
    combined_data.right_hip_std = std(all_right_hip_cycles, 0, 1);
end
if ~isempty(all_left_hip_cycles)
    combined_data.left_hip_mean = mean(all_left_hip_cycles, 1);
    combined_data.left_hip_std = std(all_left_hip_cycles, 0, 1);
end
if ~isempty(all_right_knee_cycles)
    combined_data.right_knee_mean = mean(all_right_knee_cycles, 1);
    combined_data.right_knee_std = std(all_right_knee_cycles, 0, 1);
end
if ~isempty(all_left_knee_cycles)
    combined_data.left_knee_mean = mean(all_left_knee_cycles, 1);
    combined_data.left_knee_std = std(all_left_knee_cycles, 0, 1);
end

mat_filename = sprintf('combined_gait_data_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(mat_filename, 'combined_data');

%% Print final summary
fprintf('\n=== COMBINED ANALYSIS COMPLETE ===\n');
fprintf('Files processed: %d\n', length(amc_files));
fprintf('Total cycles: %d\n', total_cycles_found);
fprintf('Right hip cycles: %d\n', size(all_right_hip_cycles, 1));
fprintf('Left hip cycles: %d\n', size(all_left_hip_cycles, 1));
fprintf('Right knee cycles: %d\n', size(all_right_knee_cycles, 1));
fprintf('Left knee cycles: %d\n', size(all_left_knee_cycles, 1));
fprintf('\nFiles saved:\n');
fprintf('  Plot: %s\n', combined_fig_filename);
fprintf('  Data: %s\n', mat_filename);

% Print cycle distribution by file
fprintf('\nCycle distribution by file:\n');
for file_idx = 1:length(amc_files)
    right_count = sum(right_hip_file_indices == file_idx);
    left_count = sum(left_hip_file_indices == file_idx);
    total_count = right_count + left_count;
    fprintf('  %s: %d cycles (R:%d, L:%d)\n', amc_files(file_idx).name, total_count, right_count, left_count);
end