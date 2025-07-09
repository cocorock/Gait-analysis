
% Process_Gait_Data_Advanced.m
close all;
clc;
clear;

% Add necessary paths
addpath('./Functions_rev/');
addpath('./Gait Data/');

% Load the all_trajectories data
load('./Gait Data/all_trajectories.mat', 'all_trajectories');

% Define parameters
sampling_freq = 120; % Hz
dt = 1 / sampling_freq; % seconds
cutoff_freq = 6; % Hz for Butterworth filter
target_resample_points = 200; % for resampling

processed_gait_data = cell(size(all_trajectories));

for i = 1:length(all_trajectories)
    fprintf('Processing trajectory %d of %d...\n', i, length(all_trajectories));
    current_trajectory_struct = all_trajectories{i};

    % 1. Detect heel strikes
    left_hs_indices = detect_heel_strikes(current_trajectory_struct.left_ankle_pos, sampling_freq);
    right_hs_indices = detect_heel_strikes(current_trajectory_struct.right_ankle_pos, sampling_freq);
    
    close all;
    figure()
    hold on
    plot(current_trajectory_struct.time, current_trajectory_struct.left_ankle_pos(:,2), 'DisplayName', 'time-Y');
    plot(current_trajectory_struct.time, current_trajectory_struct.right_ankle_pos(:,1), 'DisplayName', 'time-X');
    stem(current_trajectory_struct.time(left_hs_indices), current_trajectory_struct.left_ankle_pos(left_hs_indices,2), 'DisplayName', 'min time-Y');
    stem(current_trajectory_struct.time(right_hs_indices), current_trajectory_struct.right_ankle_pos(right_hs_indices,1), 'DisplayName', 'min time-X');
    plot(current_trajectory_struct.left_ankle_pos(:,1), current_trajectory_struct.left_ankle_pos(:,2), 'DisplayName', 'X-Y');
    stem(current_trajectory_struct.left_ankle_pos(left_hs_indices,1), current_trajectory_struct.left_ankle_pos(left_hs_indices,2), 'DisplayName', 'min X-Y');

    %FR2
    plot(current_trajectory_struct.left_ankle_pos_FR2(:,1), current_trajectory_struct.left_ankle_pos_FR2(:,2), 'DisplayName', 'FR2 X-Y');
    stem(current_trajectory_struct.left_ankle_pos_FR2(left_hs_indices,1), current_trajectory_struct.left_ankle_pos_FR2(left_hs_indices,2), 'DisplayName', 'min FR2 X-Y');
    plot(current_trajectory_struct.time, current_trajectory_struct.left_ankle_pos_FR2(:,2), 'DisplayName', 'FR2 time-Y');
    stem(current_trajectory_struct.time(left_hs_indices,1), current_trajectory_struct.left_ankle_pos_FR2(left_hs_indices,2), 'DisplayName', 'min FR2 time-Y');
    grid on;
    legend;
    
    % 2. Segment gait cycles
    % This function will return a cell array of structs, each representing a gait cycle
    segmented_cycles_for_current_trajectory = segment_gait_cycles(current_trajectory_struct, left_hs_indices, right_hs_indices);

    processed_cycles = cell(size(segmented_cycles_for_current_trajectory));

    for j = 1:length(segmented_cycles_for_current_trajectory)
        current_cycle_struct = segmented_cycles_for_current_trajectory{j};
        processed_cycle_struct = struct();
        
        field_names = fieldnames(current_cycle_struct);
        
        for k = 1:length(field_names)
            field = field_names{k};
            original_data = current_cycle_struct.(field);
            
            % Apply filtering to all fields
            filtered_data = apply_butterworth_filter(original_data, sampling_freq, cutoff_freq);
            
            % Calculate velocity for specific fields
            if contains(field, 'ankle_pos') || contains(field, 'ankle_pos_FR2')
                velocity_data = calculate_velocity(filtered_data, dt);
                processed_cycle_struct.([field '_velocity']) = velocity_data; % Store velocity
            end
            
            % Resample all fields
            % Ensure 'time' field is used for resampling, and it exists
            if isfield(current_cycle_struct, 'time')
                resampled_data = resample_trajectory(filtered_data, current_cycle_struct.time, target_resample_points);
                processed_cycle_struct.(field) = resampled_data; % Store resampled data
            else
                warning('Time field not found in current cycle struct. Skipping resampling for field %s.', field);
                processed_cycle_struct.(field) = filtered_data; % Store filtered data if no time for resampling
            end
        end
        processed_cycles{j} = processed_cycle_struct;
    end
    processed_gait_data{i} = processed_cycles;
end

% Save the processed data
output_file = './Gait Data/processed_gait_data.mat';
save(output_file, 'processed_gait_data');

fprintf('All gait data processed and saved to %s\n', output_file);
