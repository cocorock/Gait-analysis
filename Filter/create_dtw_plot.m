close all;
clc;
clear;

% Add necessary paths
addpath('../Functions_rev/');
addpath('../Gait Data/');

% --- Parameters ---
subject = '39';
filename = sprintf('../Gait Data/all_trajectories_ALL#%s.mat', subject);
sampling_freq = 120; % Hz

% --- Load Data ---
load(filename, 'all_trajectories');

% --- Extract multiple gait cycles of different durations ---
all_gait_cycles = {};
min_cycles_to_extract = 3;
max_cycles_to_extract = 6;

for i = 1:length(all_trajectories)
    trajectory = all_trajectories{i};
    heel_strikes = detect_heel_strikes(trajectory.left_ankle_pos_FR1, sampling_freq);
    
    % Ensure at least two heel strikes to form a cycle
    if length(heel_strikes) >= 2
        for j = 1:length(heel_strikes)-1
            start_frame = heel_strikes(j);
            end_frame = heel_strikes(j+1);
            
            % Extract the gait cycle data
            current_cycle_data = trajectory.left_ankle_orientation_FR1(start_frame:end_frame);
            
            % Add to the list of gait cycles
            all_gait_cycles{end+1} = current_cycle_data;
            
            % Stop extracting if max_cycles_to_extract is reached
            if length(all_gait_cycles) >= max_cycles_to_extract
                break;
            end
        end
    end
    % Stop processing trajectories if max_cycles_to_extract is reached
    if length(all_gait_cycles) >= max_cycles_to_extract
        break;
    end
end

% Filter out cycles that are too short or too long, or if there are not enough cycles
valid_gait_cycles = {};
for k = 1:length(all_gait_cycles)
    if length(all_gait_cycles{k}) > 10 % Arbitrary minimum length to avoid very short cycles
        valid_gait_cycles{end+1} = all_gait_cycles{k};
    end
end

if length(valid_gait_cycles) < min_cycles_to_extract
    error('Could not find enough valid gait cycles for DTW alignment. Found %d, need at least %d.', length(valid_gait_cycles), min_cycles_to_extract);
end

% Use only up to max_cycles_to_extract valid cycles
all_gait_cycles = valid_gait_cycles(1:min(length(valid_gait_cycles), max_cycles_to_extract));

% Choose the first cycle as the reference for DTW
reference_cycle = all_gait_cycles{1};

% Initialize arrays to store warped cycles
warped_cycles = cell(1, length(all_gait_cycles));

% Perform DTW of the reference cycle against itself to get the common warped length
[~, ix_ref, ~] = dtw(reference_cycle, reference_cycle);
common_warped_length = length(reference_cycle(ix_ref));

% Store the reference cycle, resampled to the common warped length
warped_cycles{1} = interp1((0:length(reference_cycle)-1), reference_cycle, linspace(0, length(reference_cycle)-1, common_warped_length));
warped_cycles{1} = warped_cycles{1}(:); % Ensure column vector

% --- DTW Alignment for multiple cycles ---
for k = 2:length(all_gait_cycles)
    current_cycle = all_gait_cycles{k};
    
    % Perform DTW between the reference cycle and the current cycle
    [~, ~, iy] = dtw(reference_cycle, current_cycle);
    
    % Get the warped current cycle
    warped_current_cycle = current_cycle(iy);
    
    % Resample the warped current cycle to the common warped length
    resampled_warped_current_cycle = interp1((0:length(warped_current_cycle)-1), warped_current_cycle, linspace(0, length(warped_current_cycle)-1, common_warped_length));
    
    warped_cycles{k} = resampled_warped_current_cycle(:); % Ensure column vector
end

% --- Generate Plot ---
figure('Position', [100, 100, 1200, 600]); % Adjust figure size

% Panel (a): Raw Data
subplot(1, 2, 1);
hold on;
colors = lines(length(all_gait_cycles)); % Generate distinct colors

for k = 1:length(all_gait_cycles)
    current_cycle_data = all_gait_cycles{k};
    time_vector = (0:length(current_cycle_data)-1) / sampling_freq;
    plot(time_vector, current_cycle_data, 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Gait Cycle %d', k));
end

title('(a) Raw Gait Cycles');
xlabel('Time (s)');
ylabel('Angular Position (rad)');
legend('show');
grid on;

% Panel (b): Aligned Data
subplot(1, 2, 2);
hold on;

% Determine the common normalized time base
normalized_time_base = (0:common_warped_length-1) / (common_warped_length-1);

for k = 1:length(warped_cycles)
    plot(normalized_time_base, warped_cycles{k}, 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Gait Cycle %d (Warped)', k));
end

title('(b) Aligned Gait Cycles (DTW)');
xlabel('Normalized Time');
ylabel('Angular Position (rad)');
legend('show');
grid on;

% --- Save the Plot ---
output_filename = '../Filter/dtw_normalization_plot.png';
saveas(gcf, output_filename);
fprintf('Plot saved to %s\n', output_filename);
