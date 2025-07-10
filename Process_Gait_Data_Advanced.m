

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
show_debug_plot = true;

if false
    figure(6)
    hold on;
    for i=1:size(all_trajectories)
        plot( all_trajectories{i}.left_ankle_pos(:,1) , all_trajectories{i}.left_ankle_pos(:,2));
        plot( all_trajectories{i}.right_ankle_pos(:,1) , all_trajectories{i}.right_ankle_pos(:,2)+0.1);
    end
end

for i = 1:length(all_trajectories)
    fprintf('Processing trajectory %d of %d...\n', i, length(all_trajectories));
    current_trajectory_struct = all_trajectories{i};

    % 1. Detect heel strikes
    left_hs_indices = detect_heel_strikes(current_trajectory_struct.left_ankle_pos_FR2, sampling_freq);
    right_hs_indices = detect_heel_strikes(current_trajectory_struct.right_ankle_pos_FR2, sampling_freq);
    
     if show_debug_plot
        figure(7);
        hold on
        stem(current_trajectory_struct.time(right_hs_indices), current_trajectory_struct.right_ankle_pos(right_hs_indices,1), 'DisplayName', 'min time-X');
        plot(current_trajectory_struct.left_ankle_pos(:,1), current_trajectory_struct.left_ankle_pos(:,2), '--', 'DisplayName', 'L X-Y');
        stem(current_trajectory_struct.left_ankle_pos(left_hs_indices,1), current_trajectory_struct.left_ankle_pos(left_hs_indices,2), '--', 'LineWidth', 2, 'DisplayName', 'L min X-Y');
        plot(current_trajectory_struct.right_ankle_pos(:,1), current_trajectory_struct.right_ankle_pos(:,2)+0.3, '--' ,'DisplayName', 'L X-Y');
        stem(current_trajectory_struct.right_ankle_pos(right_hs_indices,1), current_trajectory_struct.right_ankle_pos(right_hs_indices,2)+0.3, '--', 'LineWidth', 2, 'DisplayName', 'L min X-Y');
        grid on;
        legend;
    end
    
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
            
            % Calculate velocity for specific fields, then resample it
            if contains(field, 'ankle_pos') || contains(field, 'ankle_pos_FR2')
                velocity_data = calculate_velocity(filtered_data, dt);
                
                if isfield(current_cycle_struct, 'time')
                    % Create a time vector for the velocity data, which can be one point shorter
                    time_for_velocity = current_cycle_struct.time(1:size(velocity_data, 1));
                    resampled_velocity_data = resample_trajectory(velocity_data, time_for_velocity, target_resample_points);
                    processed_cycle_struct.([field '_velocity']) = resampled_velocity_data; % Store resampled velocity
                else
                    warning('Time field not found. Storing non-resampled velocity for field %s.', field);
                    processed_cycle_struct.([field '_velocity']) = velocity_data; % Store non-resampled velocity
                end
            end
            
            % Resample all fields
            % Ensure 'time' field is used for resampling, and it exists
            if isfield(current_cycle_struct, 'time')
                resampled_data = resample_trajectory(filtered_data, current_cycle_struct.time, target_resample_points);
                
                % If the current field is 'time', normalize it from 0 to 1
                if strcmp(field, 'time')
                    resampled_data = (resampled_data - min(resampled_data)) / (max(resampled_data) - min(resampled_data));
                end
                
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

processed_gait_data = {processed_gait_data{1}{:} , processed_gait_data{2}{:}};
for i=1:size(processed_gait_data,2)
    processed_gait_data{i}.time = (0:199)/200;
end
% Save the processed data
output_file = './Gait Data/new_processed_gait_data.mat';
save(output_file, 'processed_gait_data');

fprintf('All gait data processed and saved to %s\n', output_file);

% --- New Plotting Logic ---
% Flag to enable/disable these plots
show_new_plots = true;

if show_new_plots
    % Figure 1: Ankle Positions (Original Frame)
    figure('Name', 'All Ankle Positions (Original Frame)', 'Position', [100, 100, 800, 600]);
    hold on;
    title('All Ankle Positions (Original Frame)');
    xlabel('X Position');
    ylabel('Y Position');
    grid on;
%     axis equal; % Maintain aspect ratio

    % Figure 2: Ankle Velocities (Original Frame)
    figure('Name', 'All Ankle Velocities (Original Frame)', 'Position', [950, 100, 800, 600]);
    hold on;
    title('All Ankle Velocities (Original Frame)');
    xlabel('X Velocity');
    ylabel('Y Velocity');
    grid on;
%     axis equal; % Maintain aspect ratio

    % Figure 3: Ankle Positions (FR2 Frame)
    figure('Name', 'All Ankle Positions (FR2 Frame)', 'Position', [100, 750, 800, 600]);
    hold on;
    title('All Ankle Positions (FR2 Frame)');
    xlabel('X Position (FR2)');
    ylabel('Y Position (FR2)');
    grid on;
    axis equal; % Maintain aspect ratio

    % Figure 4: Ankle Velocities (FR2 Frame)
    figure('Name', 'All Ankle Velocities (FR2 Frame)', 'Position', [950, 750, 800, 600]);
    hold on;
    title('All Ankle Velocities (FR2 Frame)');
    xlabel('X Velocity (FR2)');
    ylabel('Y Velocity (FR2)');
    grid on;
    axis equal; % Maintain aspect ratio
    
    offs =0;
    % Loop through all processed gait data to plot
    for i = 1:length(processed_gait_data) % Loop through each original trajectory
%         current_trajectory_cycles = processed_gait_data{i};
%         for j = 1:length(current_trajectory_cycles) % Loop through each gait cycle within the trajectory
            current_cycle = processed_gait_data{i};

            % Plot Ankle Positions (Original Frame)
            figure(1);
%             offs = offs + 0.2;
            plot(current_cycle.ankle_pos(:,1), current_cycle.ankle_pos(:,2)+offs, 'b-', 'DisplayName', sprintf('Ankle Pos (Cycle %d-%d)', i, j));
%             plot(current_cycle.ankle_pos(:,1), current_cycle.ankle_pos(:,2)+offs, 'r-', 'DisplayName', sprintf('Right Ankle Pos (Cycle %d-%d)', i, j));
            
            % Add arrows for orientation on position plots
            arrow_skip = 10; % Plot arrow every 10 points to avoid clutter
            arrow_length = 0.03; % Adjust this for visual clarity

            if isfield(current_cycle, 'ankle_orientation') && isfield(current_cycle, 'ankle_pos')
                x = current_cycle.ankle_pos(1:arrow_skip:end, 1);
                y = current_cycle.ankle_pos(1:arrow_skip:end, 2);
                theta = current_cycle.ankle_orientation(1:arrow_skip:end);
                u = cos(theta) * arrow_length;
                v = sin(theta) * arrow_length;
                quiver(x, y+offs, u, v, 0, 'r', 'HandleVisibility', 'off');
            end

            % Plot Ankle Velocities (Original Frame)
            figure(2);
            if isfield(current_cycle, 'ankle_pos_velocity')
                plot(current_cycle.ankle_pos_velocity(:,1), current_cycle.ankle_pos_velocity(:,2), 'b-', 'DisplayName', sprintf('Left Ankle Vel (Cycle %d-%d)', i, j));
            end
            if isfield(current_cycle, 'ankle_pos_velocity')
                plot(current_cycle.ankle_pos_velocity(:,1), current_cycle.ankle_pos_velocity(:,2), 'r-', 'DisplayName', sprintf('Right Ankle Vel (Cycle %d-%d)', i, j));
            end

            % Plot Ankle Positions (FR2 Frame)
            figure(3);
            plot(current_cycle.ankle_pos_FR2(:,1), current_cycle.ankle_pos_FR2(:,2), 'b--', 'DisplayName', sprintf('Left Ankle Pos FR2 (Cycle %d-%d)', i, j));
%             plot(current_cycle.ankle_pos_FR2(:,1), current_cycle.ankle_pos_FR2(:,2), 'r--', 'DisplayName', sprintf('Right Ankle Pos FR2 (Cycle %d-%d)', i, j));
            axis equal;
            
            if isfield(current_cycle, 'ankle_orientation_FR2') && isfield(current_cycle, 'ankle_pos_FR2')
                x = current_cycle.ankle_pos_FR2(1:arrow_skip:end, 1);
                y = current_cycle.ankle_pos_FR2(1:arrow_skip:end, 2);
                theta = current_cycle.ankle_orientation_FR2(1:arrow_skip:end);
                u = cos(theta) * arrow_length;
                v = sin(theta) * arrow_length;
                quiver(x, y, u, v, 0, 'r', 'HandleVisibility', 'off');
            end

            % Plot Ankle Velocities (FR2 Frame)
            figure(4);
            if isfield(current_cycle, 'ankle_pos_FR2_velocity')
                plot(current_cycle.ankle_pos_FR2_velocity(:,1), current_cycle.ankle_pos_FR2_velocity(:,2), 'b--', 'DisplayName', sprintf('Left Ankle Vel FR2 (Cycle %d-%d)', i, j));
            end

            axis equal;
        end
    
    
    % Add legends to all figures   legend('Location', 'best'); 
    figure(1); hold off; legend('Location', 'best'); 
    figure(2); hold off; legend('Location', 'best'); 
    figure(3); hold off; legend('Location', 'best'); 
    figure(4); hold off; legend('Location', 'best'); 
   end



