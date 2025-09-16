% create_dtw_velocity_plot.m

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
all_gait_cycles_pos = {}; % Store position data first
min_cycles_to_extract = 3;
max_cycles_to_extract = 7;

for i = 1:length(all_trajectories)
    trajectory = all_trajectories{i};
    heel_strikes = detect_heel_strikes(trajectory.left_ankle_pos_FR1, sampling_freq);
    
    % Ensure at least two heel strikes to form a cycle
    if length(heel_strikes) >= 2
        for j = 1:length(heel_strikes)-1
            start_frame = heel_strikes(j);
            end_frame = heel_strikes(j+1);
            
            % Extract the gait cycle position data
            current_cycle_pos_data = trajectory.left_ankle_orientation_FR1(start_frame:end_frame);
            
            % Add to the list of gait cycles
            all_gait_cycles_pos{end+1} = current_cycle_pos_data;
            
            % Stop extracting if max_cycles_to_extract is reached
            if length(all_gait_cycles_pos) >= max_cycles_to_extract
                break;
            end
        end
    end
    % Stop processing trajectories if max_cycles_to_extract is reached
    if length(all_gait_cycles_pos) >= max_cycles_to_extract
        break;
    end
end

% Filter out cycles that are too short or too long, or if there are not enough cycles
valid_gait_cycles_pos = {};
for k = 1:length(all_gait_cycles_pos)
    if length(all_gait_cycles_pos{k}) > 10 % Arbitrary minimum length to avoid very short cycles
        valid_gait_cycles_pos{end+1} = all_gait_cycles_pos{k};
    end
end

if length(valid_gait_cycles_pos) < min_cycles_to_extract
    error('Could not find enough valid gait cycles for DTW alignment. Found %d, need at least %d.', length(valid_gait_cycles_pos), min_cycles_to_extract);
end

% Use only up to max_cycles_to_extract valid cycles
all_gait_cycles_pos = valid_gait_cycles_pos(1:min(length(valid_gait_cycles_pos), max_cycles_to_extract));

% --- Calculate Angular Velocity for each gait cycle ---
all_gait_cycles_vel = cell(1, length(all_gait_cycles_pos));
for k = 1:length(all_gait_cycles_pos)
    % Calculate angular velocity from position data
    % As seen in create_filtering_plot.m, velocity is diff(data) * sampling_freq
    current_vel = diff(all_gait_cycles_pos{k}) * sampling_freq;
    all_gait_cycles_vel{k} = current_vel(:); % Ensure column vector
end

% Choose the first velocity cycle as the reference for DTW
reference_cycle_vel = all_gait_cycles_vel{1};

% Initialize arrays to store warped velocity cycles
warped_cycles_vel = cell(1, length(all_gait_cycles_vel));

% Perform DTW of the reference cycle against itself to get the common warped length
[~, ix_ref, ~] = dtw(reference_cycle_vel, reference_cycle_vel);
common_warped_length = length(reference_cycle_vel(ix_ref));

% Store the reference cycle, resampled to the common warped length
warped_cycles_vel{1} = interp1((0:length(reference_cycle_vel)-1), reference_cycle_vel, linspace(0, length(reference_cycle_vel)-1, common_warped_length));
warped_cycles_vel{1} = warped_cycles_vel{1}(:); % Ensure column vector

% --- DTW Alignment for multiple velocity cycles ---
for k = 2:length(all_gait_cycles_vel)
    current_cycle_vel = all_gait_cycles_vel{k};
    
    % Perform DTW between the reference velocity cycle and the current velocity cycle
    [~, ~, iy] = dtw(reference_cycle_vel, current_cycle_vel);
    
    % Get the warped current velocity cycle
    warped_current_cycle_vel = current_cycle_vel(iy);
    
    % Resample the warped current velocity cycle to the common warped length
    resampled_warped_current_cycle_vel = interp1((0:length(warped_current_cycle_vel)-1), warped_current_cycle_vel, linspace(0, length(warped_current_cycle_vel)-1, common_warped_length));
    
    warped_cycles_vel{k} = resampled_warped_current_cycle_vel(:); % Ensure column vector
end

% --- Generate Plot ---
figure('Position', [100, 100, 1200, 600]); % Adjust figure size

% Panel (a): Raw Velocity Data
subplot(1, 2, 1);
hold on;
colors = lines(length(all_gait_cycles_vel)); % Generate distinct colors

for k = 1:length(all_gait_cycles_vel)
    current_cycle_data = all_gait_cycles_vel{k};
    % Velocity data is one sample shorter than position data, so adjust time vector
    time_vector = (0:length(current_cycle_data)-1) / sampling_freq;
    plot(time_vector, current_cycle_data, 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Gait Cycle %d', k));
end

title('(a) Raw Gait Cycles (Angular Velocity)');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('show');
grid on;

% Panel (b): Aligned Velocity Data
subplot(1, 2, 2);
hold on;

% Determine the common normalized time base
normalized_time_base = (0:common_warped_length-1) / (common_warped_length-1);

for k = 1:length(warped_cycles_vel)
    plot(normalized_time_base, warped_cycles_vel{k}, 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Gait Cycle %d (Warped)', k));
end

title('(b) Aligned Gait Cycles (DTW) (Angular Velocity)');
xlabel('Normalized Time');
ylabel('Angular Velocity (rad/s)');
legend('show');
grid on;

% --- Save the Plot ---
output_filename = '../Filter/dtw_normalization_velocity_plot.png';
saveas(gcf, output_filename);
fprintf('Plot saved to %s
', output_filename);
