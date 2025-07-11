
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
        plot( all_trajectories{i}.left_ankle_pos_FR1(:,1) , all_trajectories{i}.left_ankle_pos_FR1(:,2));
        plot( all_trajectories{i}.right_ankle_pos_FR1(:,1) , all_trajectories{i}.right_ankle_pos_FR1(:,2)+0.1);
    end
end

for i = 1:length(all_trajectories)
    fprintf('Processing trajectory %d of %d...\n', i, length(all_trajectories));
    current_trajectory_struct = all_trajectories{i};

    % 1. Detect heel strikes
    left_hs_indices = detect_heel_strikes(current_trajectory_struct.left_ankle_pos_FR1, sampling_freq);
    right_hs_indices = detect_heel_strikes(current_trajectory_struct.right_ankle_pos_FR1, sampling_freq);
    
     if show_debug_plot
        figure(17);
        hold on
%         stem(current_trajectory_struct.time(right_hs_indices), current_trajectory_struct.right_ankle_pos_FR1(right_hs_indices,1), 'DisplayName', 'min time-X');
        plot(current_trajectory_struct.left_ankle_pos_FR2(:,1), current_trajectory_struct.left_ankle_pos_FR2(:,2), '--', 'DisplayName', 'L X-Y');
        stem(current_trajectory_struct.left_ankle_pos_FR2(left_hs_indices,1), current_trajectory_struct.left_ankle_pos_FR2(left_hs_indices,2), '--', 'LineWidth', 2, 'DisplayName', 'L min X-Y');
        r = rand(1);
        plot(current_trajectory_struct.right_ankle_pos_FR2(:,1), current_trajectory_struct.right_ankle_pos_FR2(:,2)+r, '--' ,'DisplayName', 'L X-Y');
        stem(current_trajectory_struct.right_ankle_pos_FR2(right_hs_indices,1), current_trajectory_struct.right_ankle_pos_FR2(right_hs_indices,2)+r, '--', 'LineWidth', 2, 'DisplayName', 'L min X-Y');
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

processed_gait_data = [processed_gait_data{1}(:)' , processed_gait_data{2}(:)'];


for i=1:size(processed_gait_data, 2)
    %Normalize time
    processed_gait_data{i}.time = ((0:199)/200)';
    % reshapre rotation Matrix A
    processed_gait_data{i}.ankle_A_FR1 = reshape(processed_gait_data{i}.ankle_A_FR1, [200, 2, 2]);
    processed_gait_data{i}.ankle_A_FR2 = reshape(processed_gait_data{i}.ankle_A_FR2, [200, 2, 2]);
end

% processed_gait_data = set_lastPoint_FR1(processed_gait_data);

% % Save the processed data
output_file = './Gait Data/new_processed_gait_data.mat';
save(output_file, 'processed_gait_data');

fprintf('All gait data processed and saved to %s\n', output_file);
%  fprintf('DATA NOT saved \n');
 
% --- New Plotting Logic ---
% Flag to enable/disable these plots
show_new_plots = true;

%%

if show_new_plots
    % Create figure for phase space plots
    figure('Name', 'Phase Space Plots', 'Position', [100, 100, 1600, 1200]);

    % Subplot 1: Ankle Positions (Original Frame)
    subplot(2, 2, 1);
    hold on;
    plot(0, 0, 'k*', 'MarkerSize', 10);
    title('All Ankle Positions (Original Frame)');
    xlabel('X Position');
    ylabel('Y Position');
    grid on;
    axis equal;

    % Subplot 2: Ankle Velocities (Original Frame)
    subplot(2, 2, 2);
    hold on;
    plot(0, 0, 'k*', 'MarkerSize', 10);
    title('All Ankle Velocities (Original Frame)');
    xlabel('X Velocity');
    ylabel('Y Velocity');
    grid on;
%     axis equal;

    % Subplot 3: Ankle Positions (FR2 Frame)
    subplot(2, 2, 3);
    hold on;
    plot(0, 0, 'k*', 'MarkerSize', 10);
    title('All Ankle Positions (FR2 Frame)');
    xlabel('X Position (FR2)');
    ylabel('Y Position (FR2)');
    grid on;
    axis equal;

    % Subplot 4: Ankle Velocities (FR2 Frame)
    subplot(2, 2, 4);
    hold on;
    plot(0, 0, 'k*', 'MarkerSize', 10);
    title('All Ankle Velocities (FR2 Frame)');
    xlabel('X Velocity (FR2)');
    ylabel('Y Velocity (FR2)');
    grid on;
%     axis equal;

    % Create figure for time series plots
    figure('Name', 'Time Series Plots', 'Position', [100, 100, 1600, 1200]);

    % Subplot 5: Ankle Positions vs. Time (Original Frame)
    subplot(2, 2, 1);
    hold on;
    title('Ankle Positions vs. Time (Original Frame)');
    xlabel('Normalized Time');
    ylabel('Position');
    grid on;

    % Subplot 6: Ankle Velocities vs. Time (Original Frame)
    subplot(2, 2, 2);
    hold on;
    title('Ankle Velocities vs. Time (Original Frame)');
    xlabel('Normalized Time');
    ylabel('Velocity');
    grid on;

    % Subplot 7: Ankle Positions vs. Time (FR2 Frame)
    subplot(2, 2, 3);
    hold on;
    title('Ankle Positions vs. Time (FR2 Frame)');
    xlabel('Normalized Time');
    ylabel('Position (FR2)');
    grid on;

    % Subplot 8: Ankle Velocities vs. Time (FR2 Frame)
    subplot(2, 2, 4);
    hold on;
    title('Ankle Velocities vs. Time (FR2 Frame)');
    xlabel('Normalized Time');
    ylabel('Velocity (FR2)');
    grid on;
    
    
    offs =0;
    % Loop through all processed gait data to plot
    for i = 1:length(processed_gait_data) % Loop through each original trajectory
%         current_trajectory_cycles = processed_gait_data{i};
%         for j = 1:length(current_trajectory_cycles) % Loop through each gait cycle within the trajectory
            current_cycle = processed_gait_data{i};

            % --- Phase Space Plots ---
            figure(1);

            % Plot Ankle Positions (Original Frame)
            subplot(2, 2, 1);
            plot(current_cycle.ankle_pos_FR1(:,1), current_cycle.ankle_pos_FR1(:,2)+offs, 'b-', 'DisplayName', sprintf('Ankle Pos (Cycle %d)', i));
            plot(current_cycle.ankle_pos_FR1(1,1), current_cycle.ankle_pos_FR1(1,2)+offs, 'r+', 'MarkerSize', 10, 'HandleVisibility', 'off');
            if isfield(current_cycle, 'ankle_orientation') && isfield(current_cycle, 'ankle_pos')
                arrow_skip = 10;
                arrow_length = 0.03;
                x = current_cycle.ankle_pos_FR1(1:arrow_skip:end, 1);
                y = current_cycle.ankle_pos_FR1(1:arrow_skip:end, 2);
                theta = current_cycle.ankle_orientation(1:arrow_skip:end);
                u = cos(theta) * arrow_length;
                v = sin(theta) * arrow_length;
                quiver(x, y+offs, u, v, 0, 'r', 'HandleVisibility', 'off');
            end

            % Plot Ankle Velocities (Original Frame)
            subplot(2, 2, 2);
            if isfield(current_cycle, 'ankle_pos_FR1_velocity')
                plot(current_cycle.ankle_pos_FR1_velocity(:,1), current_cycle.ankle_pos_FR1_velocity(:,2), 'b-', 'DisplayName', sprintf('Left Ankle Vel (Cycle %d)', i));
                plot(current_cycle.ankle_pos_FR1_velocity(1,1), current_cycle.ankle_pos_FR1_velocity(1,2), 'r+', 'MarkerSize', 10, 'HandleVisibility', 'off');
            end

            % Plot Ankle Positions (FR2 Frame)
            subplot(2, 2, 3);
            plot(current_cycle.ankle_pos_FR2(:,1), current_cycle.ankle_pos_FR2(:,2), 'b--', 'DisplayName', sprintf('Left Ankle Pos FR2 (Cycle %d)', i));
            plot(current_cycle.ankle_pos_FR2(1,1), current_cycle.ankle_pos_FR2(1,2), 'r+', 'MarkerSize', 10, 'HandleVisibility', 'off');
            if isfield(current_cycle, 'ankle_orientation_FR2') && isfield(current_cycle, 'ankle_pos_FR2')
                arrow_skip = 10;
                arrow_length = 0.03;
                x = current_cycle.ankle_pos_FR2(1:arrow_skip:end, 1);
                y = current_cycle.ankle_pos_FR2(1:arrow_skip:end, 2);
                theta = current_cycle.ankle_orientation_FR2(1:arrow_skip:end);
                u = cos(theta) * arrow_length;
                v = sin(theta) * arrow_length;
                quiver(x, y, u, v, 0, 'r', 'HandleVisibility', 'off');
            end

            % Plot Ankle Velocities (FR2 Frame)
            subplot(2, 2, 4);
            if isfield(current_cycle, 'ankle_pos_FR2_velocity')
                plot(current_cycle.ankle_pos_FR2_velocity(:,1), current_cycle.ankle_pos_FR2_velocity(:,2), 'b--', 'DisplayName', sprintf('Left Ankle Vel FR2 (Cycle %d)', i));
                plot(current_cycle.ankle_pos_FR2_velocity(1,1), current_cycle.ankle_pos_FR2_velocity(1,2), 'r+', 'MarkerSize', 10, 'HandleVisibility', 'off');
            end

            % --- Time Series Plots ---
            figure(2);

            % Plot Ankle Positions vs. Time (Original Frame)
            subplot(2, 2, 1);
            plot(current_cycle.time, current_cycle.ankle_pos_FR1(:,1), 'b-', 'DisplayName', sprintf('X Pos (Cycle %d)', i));
            plot(current_cycle.time, current_cycle.ankle_pos_FR1(:,2), 'r-', 'DisplayName', sprintf('Y Pos (Cycle %d)', i));

            % Plot Ankle Velocities vs. Time (Original Frame)
            subplot(2, 2, 2);
            if isfield(current_cycle, 'ankle_pos_FR2_velocity')
                plot(current_cycle.time, current_cycle.ankle_pos_FR1_velocity(:,1), 'b-', 'DisplayName', sprintf('X Vel FR1 (Cycle %d)', i));
                plot(current_cycle.time, current_cycle.ankle_pos_FR1_velocity(:,2), 'r-', 'DisplayName', sprintf('Y Vel FR1 (Cycle %d)', i));
            end

            % Plot Ankle Positions vs. Time (FR2 Frame)
            subplot(2, 2, 3);
            plot(current_cycle.time, current_cycle.ankle_pos_FR2(:,1), 'b--', 'DisplayName', sprintf('X Pos FR2 (Cycle %d)', i));
            plot(current_cycle.time, current_cycle.ankle_pos_FR2(:,2), 'r--', 'DisplayName', sprintf('Y Pos FR2 (Cycle %d)', i));

            % Plot Ankle Velocities vs. Time (FR2 Frame)
            subplot(2, 2, 4);
            if isfield(current_cycle, 'ankle_pos_FR2_velocity')
                plot(current_cycle.time, current_cycle.ankle_pos_FR2_velocity(:,1), 'b--', 'DisplayName', sprintf('X Vel FR2 (Cycle %d)', i));
                plot(current_cycle.time, current_cycle.ankle_pos_FR2_velocity(:,2), 'r--', 'DisplayName', sprintf('Y Vel FR2 (Cycle %d)', i));
            end
    end
    
    
    % Add legends to all subplots
    figure(1);
    for k = 1:4
        subplot(2, 2, k);
        hold off;
        legend('Location', 'best');
    end

    figure(2);
    for k = 1:4
        subplot(2, 2, k);
        hold off;
        legend('Location', 'best');
    end
   end



