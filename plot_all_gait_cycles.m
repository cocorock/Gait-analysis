% Plot all gait cycles together using robust knee minima detection
% This script extracts gait cycles and plots them overlaid for comparison

%% Extract gait cycles using robust method
fprintf('Extracting gait cycles using robust knee minima detection...\n');
gait_cycles_data = extract_gait_cycles_knee_minima_robust('39_01.amc');

if isempty(gait_cycles_data)
    error('No gait cycles detected! Check the detection parameters.');
end

fprintf('Found %d validated gait cycles\n', length(gait_cycles_data));

%% Separate cycles by leg
right_cycles = find(strcmp({gait_cycles_data.leg}, 'right'));
left_cycles = find(strcmp({gait_cycles_data.leg}, 'left'));

fprintf('Right leg cycles: %d\n', length(right_cycles));
fprintf('Left leg cycles: %d\n', length(left_cycles));

%% Create comprehensive plot
figure('Name', 'All Gait Cycles Overlaid', 'Position', [100, 100, 1600, 1200]);

% Define colors for cycles
colors_right = lines(length(right_cycles));
colors_left = lines(length(left_cycles));

%% Plot 1: Right Hip Flexion - All Cycles
subplot(2,3,1);
hold on;
for i = 1:length(right_cycles)
    cycle_idx = right_cycles(i);
    plot(gait_cycles_data(cycle_idx).time_normalized, ...
         gait_cycles_data(cycle_idx).right_hip_flex, ...
         'Color', colors_right(i,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('R-Cycle %d', gait_cycles_data(cycle_idx).cycle_number));
end
xlabel('Normalized Time (0-1)');
ylabel('Right Hip Flexion (degrees)');
title(sprintf('Right Hip Flexion - All %d Cycles', length(right_cycles)));
grid on;
legend('Location', 'best', 'FontSize', 8);

%% Plot 2: Left Hip Flexion - All Cycles  
subplot(2,3,2);
hold on;
for i = 1:length(left_cycles)
    cycle_idx = left_cycles(i);
    plot(gait_cycles_data(cycle_idx).time_normalized, ...
         gait_cycles_data(cycle_idx).left_hip_flex, ...
         'Color', colors_left(i,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('L-Cycle %d', gait_cycles_data(cycle_idx).cycle_number));
end
xlabel('Normalized Time (0-1)');
ylabel('Left Hip Flexion (degrees)');
title(sprintf('Left Hip Flexion - All %d Cycles', length(left_cycles)));
grid on;
legend('Location', 'best', 'FontSize', 8);

%% Plot 3: Both Hip Flexion - Mean ± STD
subplot(2,3,3);
hold on;

% Calculate mean and std for right hip (from right leg cycles)
if ~isempty(right_cycles)
    % Interpolate all right cycles to same length for averaging
    interp_length = 101; % 0 to 100% in 1% steps
    right_hip_matrix = zeros(length(right_cycles), interp_length);

    for i = 1:length(right_cycles)
        cycle_idx = right_cycles(i);
        time_norm = gait_cycles_data(cycle_idx).time_normalized;
        hip_data = gait_cycles_data(cycle_idx).right_hip_flex;

        % Interpolate to standard length
        time_standard = linspace(0, 1, interp_length);
        right_hip_matrix(i, :) = interp1(time_norm, hip_data, time_standard, 'linear');
    end

    right_hip_mean = mean(right_hip_matrix, 1);
    right_hip_std = std(right_hip_matrix, 0, 1);

    time_standard = linspace(0, 1, interp_length);

    % Plot mean ± std
    fill([time_standard, fliplr(time_standard)], ...
         [right_hip_mean + right_hip_std, fliplr(right_hip_mean - right_hip_std)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Right ± STD');
    plot(time_standard, right_hip_mean, 'r-', 'LineWidth', 3, 'DisplayName', 'Right Mean');
end

% Calculate mean and std for left hip (from left leg cycles)
if ~isempty(left_cycles)
    left_hip_matrix = zeros(length(left_cycles), interp_length);

    for i = 1:length(left_cycles)
        cycle_idx = left_cycles(i);
        time_norm = gait_cycles_data(cycle_idx).time_normalized;
        hip_data = gait_cycles_data(cycle_idx).left_hip_flex;

        time_standard = linspace(0, 1, interp_length);
        left_hip_matrix(i, :) = interp1(time_norm, hip_data, time_standard, 'linear');
    end

    left_hip_mean = mean(left_hip_matrix, 1);
    left_hip_std = std(left_hip_matrix, 0, 1);

    % Plot mean ± std
    fill([time_standard, fliplr(time_standard)], ...
         [left_hip_mean + left_hip_std, fliplr(left_hip_mean - left_hip_std)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Left ± STD');
    plot(time_standard, left_hip_mean, 'b-', 'LineWidth', 3, 'DisplayName', 'Left Mean');
end

xlabel('Normalized Time (0-1)');
ylabel('Hip Flexion (degrees)');
title('Hip Flexion - Mean ± STD');
legend('Location', 'best');
grid on;

%% Plot 4: Right Knee Flexion - All Cycles
subplot(2,3,4);
hold on;
for i = 1:length(right_cycles)
    cycle_idx = right_cycles(i);
    plot(gait_cycles_data(cycle_idx).time_normalized, ...
         gait_cycles_data(cycle_idx).right_knee_flex, ...
         'Color', colors_right(i,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('R-Cycle %d', gait_cycles_data(cycle_idx).cycle_number));
end
xlabel('Normalized Time (0-1)');
ylabel('Right Knee Flexion (degrees)');
title(sprintf('Right Knee Flexion - All %d Cycles', length(right_cycles)));
grid on;
legend('Location', 'best', 'FontSize', 8);

%% Plot 5: Left Knee Flexion - All Cycles
subplot(2,3,5);
hold on;
for i = 1:length(left_cycles)
    cycle_idx = left_cycles(i);
    plot(gait_cycles_data(cycle_idx).time_normalized, ...
         gait_cycles_data(cycle_idx).left_knee_flex, ...
         'Color', colors_left(i,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('L-Cycle %d', gait_cycles_data(cycle_idx).cycle_number));
end
xlabel('Normalized Time (0-1)');
ylabel('Left Knee Flexion (degrees)');
title(sprintf('Left Knee Flexion - All %d Cycles', length(left_cycles)));
grid on;
legend('Location', 'best', 'FontSize', 8);

%% Plot 6: Both Knee Flexion - Mean ± STD
subplot(2,3,6);
hold on;

% Calculate mean and std for right knee (from right leg cycles)
if ~isempty(right_cycles)
    right_knee_matrix = zeros(length(right_cycles), interp_length);

    for i = 1:length(right_cycles)
        cycle_idx = right_cycles(i);
        time_norm = gait_cycles_data(cycle_idx).time_normalized;
        knee_data = gait_cycles_data(cycle_idx).right_knee_flex;

        time_standard = linspace(0, 1, interp_length);
        right_knee_matrix(i, :) = interp1(time_norm, knee_data, time_standard, 'linear');
    end

    right_knee_mean = mean(right_knee_matrix, 1);
    right_knee_std = std(right_knee_matrix, 0, 1);

    % Plot mean ± std
    fill([time_standard, fliplr(time_standard)], ...
         [right_knee_mean + right_knee_std, fliplr(right_knee_mean - right_knee_std)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Right ± STD');
    plot(time_standard, right_knee_mean, 'r-', 'LineWidth', 3, 'DisplayName', 'Right Mean');
end

% Calculate mean and std for left knee (from left leg cycles)
if ~isempty(left_cycles)
    left_knee_matrix = zeros(length(left_cycles), interp_length);

    for i = 1:length(left_cycles)
        cycle_idx = left_cycles(i);
        time_norm = gait_cycles_data(cycle_idx).time_normalized;
        knee_data = gait_cycles_data(cycle_idx).left_knee_flex;

        time_standard = linspace(0, 1, interp_length);
        left_knee_matrix(i, :) = interp1(time_norm, knee_data, time_standard, 'linear');
    end

    left_knee_mean = mean(left_knee_matrix, 1);
    left_knee_std = std(left_knee_matrix, 0, 1);

    % Plot mean ± std
    fill([time_standard, fliplr(time_standard)], ...
         [left_knee_mean + left_knee_std, fliplr(left_knee_mean - left_knee_std)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Left ± STD');
    plot(time_standard, left_knee_mean, 'b-', 'LineWidth', 3, 'DisplayName', 'Left Mean');
end

xlabel('Normalized Time (0-1)');
ylabel('Knee Flexion (degrees)');
title('Knee Flexion - Mean ± STD');
legend('Location', 'best');
grid on;

%% Print summary statistics
fprintf('\n=== GAIT CYCLE ANALYSIS SUMMARY ===\n');

if ~isempty(right_cycles)
    fprintf('\nRIGHT LEG CYCLES (%d total):\n', length(right_cycles));
    right_durations = [gait_cycles_data(right_cycles).duration];
    fprintf('  Duration: %.3f ± %.3f s (range: %.3f-%.3f s)\n', ...
            mean(right_durations), std(right_durations), ...
            min(right_durations), max(right_durations));

    fprintf('  Right Hip Range: %.1f ± %.1f° (%.1f° to %.1f°)\n', ...
            mean(max(right_hip_matrix, [], 2) - min(right_hip_matrix, [], 2)), ...
            std(max(right_hip_matrix, [], 2) - min(right_hip_matrix, [], 2)), ...
            min(right_hip_matrix(:)), max(right_hip_matrix(:)));

    fprintf('  Right Knee Range: %.1f ± %.1f° (%.1f° to %.1f°)\n', ...
            mean(max(right_knee_matrix, [], 2) - min(right_knee_matrix, [], 2)), ...
            std(max(right_knee_matrix, [], 2) - min(right_knee_matrix, [], 2)), ...
            min(right_knee_matrix(:)), max(right_knee_matrix(:)));
end

if ~isempty(left_cycles)
    fprintf('\nLEFT LEG CYCLES (%d total):\n', length(left_cycles));
    left_durations = [gait_cycles_data(left_cycles).duration];
    fprintf('  Duration: %.3f ± %.3f s (range: %.3f-%.3f s)\n', ...
            mean(left_durations), std(left_durations), ...
            min(left_durations), max(left_durations));

    fprintf('  Left Hip Range: %.1f ± %.1f° (%.1f° to %.1f°)\n', ...
            mean(max(left_hip_matrix, [], 2) - min(left_hip_matrix, [], 2)), ...
            std(max(left_hip_matrix, [], 2) - min(left_hip_matrix, [], 2)), ...
            min(left_hip_matrix(:)), max(left_hip_matrix(:)));

    fprintf('  Left Knee Range: %.1f ± %.1f° (%.1f° to %.1f°)\n', ...
            mean(max(left_knee_matrix, [], 2) - min(left_knee_matrix, [], 2)), ...
            std(max(left_knee_matrix, [], 2) - min(left_knee_matrix, [], 2)), ...
            min(left_knee_matrix(:)), max(left_knee_matrix(:)));
end

fprintf('\nAll cycles normalized to 0-1 time and overlaid for comparison.\n');
fprintf('Mean ± STD plots show average gait pattern with variability.\n');