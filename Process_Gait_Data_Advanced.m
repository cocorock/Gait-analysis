
% Process_Gait_Data_Advanced.m
close all;
clc;
clear;

% Add necessary paths
addpath('./Functions_rev/');
addpath('./Gait Data/');

subject = '39'; %07 %39*
filename = sprintf('./Gait Data/all_trajectories_ALL-10#%s.mat', subject);

% Load the all_trajectories data
load(filename, 'all_trajectories');

% Define parameters
sampling_freq = 120; % Hz
dt = 1 / sampling_freq; % seconds
cutoff_freq = 6; % Hz for Butterworth filter
target_resample_points = 200; % for resampling

p1_processed_gait_data = cell(size(all_trajectories));
show_debug_plot = true;

if false
    figure(6)
    hold on;
    for i=1:size(all_trajectories)
        plot( all_trajectories{i}.left_ankle_pos_FR1(:,1) , all_trajectories{i}.left_ankle_pos_FR1(:,2));
        plot( all_trajectories{i}.right_ankle_pos_FR1(:,1) , all_trajectories{i}.right_ankle_pos_FR1(:,2)+0.1);
    end
end
voffs = 0.2;
for i = 1:length(all_trajectories)
    fprintf('\nProcessing trajectory %d of %d...\n', i, length(all_trajectories));
    current_trajectory_struct = all_trajectories{i};

    % 1. Detect heel strikes
    left_hs_indices = detect_heel_strikes(current_trajectory_struct.left_ankle_pos_FR1, sampling_freq);
    right_hs_indices = detect_heel_strikes(current_trajectory_struct.right_ankle_pos_FR1, sampling_freq);
    
     if show_debug_plot
        figure(17);
        hold on
        plot(current_trajectory_struct.left_ankle_pos_FR2(:,1), current_trajectory_struct.left_ankle_pos_FR2(:,2)+voffs, '--', 'DisplayName', 'L X-Y');
        stem(current_trajectory_struct.left_ankle_pos_FR2(left_hs_indices,1), current_trajectory_struct.left_ankle_pos_FR2(left_hs_indices,2)+voffs, '--', 'LineWidth', 2, 'DisplayName', 'L min X-Y');
        voffs = 0.2 + voffs;
        plot(current_trajectory_struct.right_ankle_pos_FR2(:,1), current_trajectory_struct.right_ankle_pos_FR2(:,2)+voffs, '--' ,'DisplayName', 'L X-Y');
        stem(current_trajectory_struct.right_ankle_pos_FR2(right_hs_indices,1), current_trajectory_struct.right_ankle_pos_FR2(right_hs_indices,2)+voffs, '--', 'LineWidth', 2, 'DisplayName', 'L min X-Y');
        voffs = 0.2 + voffs;
        grid on;
        legend;
    end
    
    % 2. Segment gait cycles
    % This function will return a cell array of structs, each representing a gait cycle
    segmented_cycles_for_current_trajectory = segment_gait_cycles(current_trajectory_struct, left_hs_indices, right_hs_indices);
    processed_cycles = cell(size(segmented_cycles_for_current_trajectory));
    
    fprintf('procesados por segment Size %d x %d \n', size(processed_cycles))
    for j = 1:length(segmented_cycles_for_current_trajectory)
        current_cycle_struct = segmented_cycles_for_current_trajectory{j};
        processed_cycle_struct = struct();

        % Subtracting refPoint from FR2 and updatating ankle_b_FR2
%         fprintf('Bef last point of ankle_pos_FR2 %2.2f, %2.2f \n', current_cycle_struct.ankle_pos_FR2(end , :));
%         fprintf('\tRefPoint %2.2f, %2.2f \n', current_cycle_struct.RefPoint);
        current_cycle_struct.ankle_pos_FR2 = current_cycle_struct.ankle_pos_FR2 - current_cycle_struct.RefPoint;
        current_cycle_struct.ankle_b_FR2 = current_cycle_struct.ankle_b_FR2 - current_cycle_struct.RefPoint;
        current_cycle_struct = rmfield(current_cycle_struct, 'RefPoint');
%         fprintf('Af last point of ankle_pos_FR2 %2.2f, %2.2f \n', current_cycle_struct.ankle_pos_FR2(end , :));

        field_names = fieldnames(current_cycle_struct);
        
        for k = 1:length(field_names)
            field = field_names{k};
            original_data = current_cycle_struct.(field);
            

            if contains(field, '_pos_') || contains(field, 'orientation')
                % The original code had a suspected typo `original_data{field}`, which is corrected to `original_data`.
                % Also, the plot flag is set to false to avoid generating excessive plots.
                filtered_data = apply_butterworth_filter(original_data, sampling_freq, cutoff_freq, false);
            else
                % For other fields (like 'time' or rotation matrices), do not apply the filter.
                filtered_data = original_data;
            end
            
            % Calculate velocity for specific fields, then resample it
            % This now uses 'filtered_data', which is either the filtered or original data.
            if contains(field, 'ankle_pos_FR1') || contains(field, 'ankle_pos_FR2')
                
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
            
            % Resample all fields using the potentially filtered data
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
    p1_processed_gait_data{i} = processed_cycles;
end

nfiles = size(p1_processed_gait_data, 1);
idx =1;
for i=1:nfiles
    gait_cycle = p1_processed_gait_data{i};
    ncycles = size(gait_cycle, 2);
    for j=1:ncycles
        processed_gait_data{idx}= gait_cycle{j};
        idx = idx + 1;
    end
    
end

fprintf('\t\t Size of processed_gait_data %d \n', length(processed_gait_data));

for i=1:size(processed_gait_data, 2)
    fprintf('\t Processing processed_gait_data %d of %d...\n', i, length(processed_gait_data));
    %Normalize time
%     processed_gait_data{i}.time = ((0:199)/200)';
    % reshapre rotation Matrix A
    processed_gait_data{i}.ankle_A_FR1 = reshape(processed_gait_data{i}.ankle_A_FR1, [200, 2, 2]);
    processed_gait_data{i}.ankle_A_FR2 = reshape(processed_gait_data{i}.ankle_A_FR2, [200, 2, 2]);
end

% processed_gait_data = set_lastPoint_FR1(processed_gait_data);

% % Save the processed data
output_file = sprintf('./Gait Data/new_processed_gait_data#%s_%d.mat', subject,size(processed_gait_data, 2));
save(output_file, 'processed_gait_data');

fprintf('All gait data processed and saved to %s\n', output_file);
%  fprintf('DATA NOT saved \n');
 
% --- New Plotting Logic ---
% Flag to enable/disable these plots
show_new_plots = true;

%%

if show_new_plots
    % Create figure for phase space plots
    figure('Name', 'Phase Space Plots', 'Position', [100, 100, 1800, 800]);

    % Row 1: Original Frame (FR1)
    subplot(2, 3, 1); hold on; title('Ankle Positions (FR1)'); xlabel('X Position'); ylabel('Y Position'); grid on; axis equal;
    subplot(2, 3, 2); hold on; title('Ankle Velocities (FR1)'); xlabel('X Velocity'); ylabel('Y Velocity'); grid on;
    subplot(2, 3, 3); hold on; title('Ankle Orientation (FR1)'); xlabel('Time (s)'); ylabel('Orientation (rad)'); grid on;

    % Row 2: Second Frame (FR2)
    subplot(2, 3, 4); hold on; title('Ankle Positions (FR2)'); xlabel('X Position (FR2)'); ylabel('Y Position (FR2)'); grid on; axis equal;
    subplot(2, 3, 5); hold on; title('Ankle Velocities (FR2)'); xlabel('X Velocity (FR2)'); ylabel('Y Velocity (FR2)'); grid on;
    subplot(2, 3, 6); hold on; title('Ankle Orientation (FR2)'); xlabel('Time (s)'); ylabel('Orientation (rad)'); grid on;

    % Create figure for time series plots
    figure('Name', 'Time Series Plots', 'Position', [150, 150, 1800, 800]);

    % Row 1: Original Frame (FR1)
    subplot(2, 3, 1); hold on; title('Ankle Positions vs. Time (FR1)'); xlabel('Normalized Time'); ylabel('Position'); grid on;
    subplot(2, 3, 2); hold on; title('Ankle Velocities vs. Time (FR1)'); xlabel('Normalized Time'); ylabel('Velocity'); grid on;
    subplot(2, 3, 3); hold on; title('Ankle Orientation vs. Time (FR1)'); xlabel('Normalized Time'); ylabel('Orientation (rad)'); grid on;

    % Row 2: Second Frame (FR2)
    subplot(2, 3, 4); hold on; title('Ankle Positions vs. Time (FR2)'); xlabel('Normalized Time'); ylabel('Position (FR2)'); grid on;
    subplot(2, 3, 5); hold on; title('Ankle Velocities vs. Time (FR2)'); xlabel('Normalized Time'); ylabel('Velocity (FR2)'); grid on;
    subplot(2, 3, 6); hold on; title('Ankle Orientation vs. Time (FR2)'); xlabel('Normalized Time'); ylabel('Orientation (rad)'); grid on;

    % Create figure for ankle bias plots
    figure('Name', 'Ankle Bias (b) and Orientation (A) Plots', 'Position', [200, 200, 1200, 1000]);
    subplot(2, 2, 1); hold on; title('Ankle Bias (b) vs. Time (FR1)'); xlabel('Normalized Time'); ylabel('Bias Value'); grid on;
    subplot(2, 2, 2); hold on; title('Ankle Bias (b) vs. Time (FR2)'); xlabel('Normalized Time'); ylabel('Bias Value'); grid on;
    subplot(2, 2, 3); hold on; title('Ankle Orientation (A) vs. Time (FR1)'); xlabel('Normalized Time'); ylabel('Angle (deg)'); grid on;
    subplot(2, 2, 4); hold on; title('Ankle Orientation (A) vs. Time (FR2)'); xlabel('Normalized Time'); ylabel('Angle (deg)'); grid on;

    % Loop through all processed gait data to plot
    for i = 1:length(processed_gait_data)
        current_cycle = processed_gait_data{i};

        % --- Phase Space Plots ---
        figure(1);
        % FR1
        subplot(2, 3, 1);
        plot(current_cycle.ankle_pos_FR1(:,1), current_cycle.ankle_pos_FR1(:,2), '--', 'HandleVisibility','off');
        plot(current_cycle.ankle_pos_FR1(1,1), current_cycle.ankle_pos_FR1(1,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));

        if isfield(current_cycle, 'ankle_pos_FR1_velocity')
            subplot(2, 3, 2);
            plot(current_cycle.ankle_pos_FR1_velocity(:,1), current_cycle.ankle_pos_FR1_velocity(:,2), '--', 'HandleVisibility','off');
            plot(current_cycle.ankle_pos_FR1_velocity(1,1), current_cycle.ankle_pos_FR1_velocity(1,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
        if isfield(current_cycle, 'ankle_orientation_FR1')
            subplot(2, 3, 3);
            plot(current_cycle.time, current_cycle.ankle_orientation_FR1, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_orientation_FR1(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
        % FR2
        subplot(2, 3, 4);
        plot(current_cycle.ankle_pos_FR2(:,1), current_cycle.ankle_pos_FR2(:,2), '--', 'HandleVisibility','off');
        plot(current_cycle.ankle_pos_FR2(1,1), current_cycle.ankle_pos_FR2(1,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));

        if isfield(current_cycle, 'ankle_pos_FR2_velocity')
            subplot(2, 3, 5);
            plot(current_cycle.ankle_pos_FR2_velocity(:,1), current_cycle.ankle_pos_FR2_velocity(:,2), '--', 'HandleVisibility','off');
            plot(current_cycle.ankle_pos_FR2_velocity(1,1), current_cycle.ankle_pos_FR2_velocity(1,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
        if isfield(current_cycle, 'ankle_orientation_FR2')
            subplot(2, 3, 6);
            plot(current_cycle.time, current_cycle.ankle_orientation_FR2, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_orientation_FR2(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        % --- Time Series Plots ---
        figure(2);
        % FR1
        subplot(2, 3, 1);
        plot(current_cycle.time, current_cycle.ankle_pos_FR1, '--', 'HandleVisibility','off');
        plot(current_cycle.time(1), current_cycle.ankle_pos_FR1(1,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));

        if isfield(current_cycle, 'ankle_pos_FR1_velocity')
            subplot(2, 3, 2);
            plot(current_cycle.time, current_cycle.ankle_pos_FR1_velocity, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_pos_FR1_velocity(1,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
        if isfield(current_cycle, 'ankle_orientation_FR1')
            subplot(2, 3, 3);
            plot(current_cycle.time, current_cycle.ankle_orientation_FR1, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_orientation_FR1(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
        % FR2
        subplot(2, 3, 4);
        plot(current_cycle.time, current_cycle.ankle_pos_FR2, '--', 'HandleVisibility','off');
        plot(current_cycle.time(1), current_cycle.ankle_pos_FR2(1,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));

        if isfield(current_cycle, 'ankle_pos_FR2_velocity')
            subplot(2, 3, 5);
            plot(current_cycle.time, current_cycle.ankle_pos_FR2_velocity, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_pos_FR2_velocity(1,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
        if isfield(current_cycle, 'ankle_orientation_FR2')
            subplot(2, 3, 6);
            plot(current_cycle.time, current_cycle.ankle_orientation_FR2, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_orientation_FR2(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        % --- Ankle Bias Plots ---
        figure(3);
        if isfield(current_cycle, 'ankle_b_FR1')
            subplot(2, 2, 1);
            stem(current_cycle.time, current_cycle.ankle_b_FR1, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_b_FR1(1,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
        if isfield(current_cycle, 'ankle_b_FR2')
            subplot(2, 2, 2);
            plot(current_cycle.time, current_cycle.ankle_b_FR2, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_b_FR2(1,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
        if isfield(current_cycle, 'ankle_A_FR1')
            subplot(2, 2, 3);
            angle_A_FR1 = atan2d(squeeze(current_cycle.ankle_A_FR1(:, 2, 1)), squeeze(current_cycle.ankle_A_FR1(:, 1, 1)));
            plot(current_cycle.time, angle_A_FR1, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), angle_A_FR1(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
        if isfield(current_cycle, 'ankle_A_FR2')
            subplot(2, 2, 4);
            angle_A_FR2 = atan2d(squeeze(current_cycle.ankle_A_FR2(:, 2, 1)), squeeze(current_cycle.ankle_A_FR2(:, 1, 1)));
            plot(current_cycle.time, angle_A_FR2, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), angle_A_FR2(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        
    end

    % Add origin marker and legends to all subplots
    for k = 1:6
        figure(1);
        subplot(2, 3, k);
        plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
%         legend('show', 'Location', 'best');
        
        figure(2);
        subplot(2, 3, k);
        plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
%         legend('show', 'Location', 'best');
    end

    % Add origin marker and legends to the new figure
    figure(3);
    for k = 1:4
        subplot(2, 2, k);
        plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
%         legend('show', 'Location', 'best');
    end

    % Create figure for transformed time series plots
    figure(4);
    set(gcf, 'Name', 'Transformed Time Series Plots', 'Position', [250, 250, 1800, 800]);

    % Row 1: Original Frame (FR1) - Transformed
    subplot(2, 3, 1); hold on; title('Transformed Ankle Positions vs. Time (FR1)'); xlabel('Normalized Time'); ylabel('Position'); grid on;
    subplot(2, 3, 2); hold on; title('Transformed Ankle Velocities vs. Time (FR1)'); xlabel('Normalized Time'); ylabel('Velocity'); grid on;
    subplot(2, 3, 3); hold on; title('Ankle Orientation vs. Time (FR1)'); xlabel('Normalized Time'); ylabel('Orientation (rad)'); grid on;

    % Row 2: Second Frame (FR2) - Transformed
    subplot(2, 3, 4); hold on; title('Transformed Ankle Positions vs. Time (FR2)'); xlabel('Normalized Time'); ylabel('Position (FR2)'); grid on;
    subplot(2, 3, 5); hold on; title('Transformed Ankle Velocities vs. Time (FR2)'); xlabel('Normalized Time'); ylabel('Velocity (FR2)'); grid on;
    subplot(2, 3, 6); hold on; title('Ankle Orientation vs. Time (FR2)'); xlabel('Normalized Time'); ylabel('Orientation (rad)'); grid on;

    % Loop through all processed gait data to plot transformed data
    for i = 1:length(processed_gait_data)
        current_cycle = processed_gait_data{i};

        % --- Transformed Time Series Plots ---
        figure(4);
        
        % FR1 Transformation and Plotting
        if isfield(current_cycle, 'ankle_pos_FR1') && isfield(current_cycle, 'ankle_A_FR1') && isfield(current_cycle, 'ankle_b_FR1')
            transformed_pos_FR1 = zeros(size(current_cycle.ankle_pos_FR1));
            for t = 1:size(current_cycle.ankle_pos_FR1, 1)
                A_inv_FR1 = squeeze(current_cycle.ankle_A_FR1(t, :, :))'; % Inverse of rotation matrix is its transpose
                b_FR1 = current_cycle.ankle_b_FR1(t, :)';
                pos_FR1 = current_cycle.ankle_pos_FR1(t, :)';
                transformed_pos_FR1(t, :) = (A_inv_FR1 * (pos_FR1 - b_FR1))';
            end
            subplot(2, 3, 1);
            plot(current_cycle.time, transformed_pos_FR1, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), transformed_pos_FR1(1,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        if isfield(current_cycle, 'ankle_pos_FR1_velocity') && isfield(current_cycle, 'ankle_A_FR1')
            transformed_vel_FR1 = zeros(size(current_cycle.ankle_pos_FR1_velocity));
            for t = 1:size(current_cycle.ankle_pos_FR1_velocity, 1)
                A_inv_FR1 = squeeze(current_cycle.ankle_A_FR1(t, :, :))';
                vel_FR1 = current_cycle.ankle_pos_FR1_velocity(t, :)';
                transformed_vel_FR1(t, :) = (A_inv_FR1 * vel_FR1)';
            end
            subplot(2, 3, 2);
            plot(current_cycle.time, transformed_vel_FR1, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), transformed_vel_FR1(1,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        if isfield(current_cycle, 'ankle_orientation_FR1')
            subplot(2, 3, 3);
            plot(current_cycle.time, current_cycle.ankle_orientation_FR1, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_orientation_FR1(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        % FR2 Transformation and Plotting
        if isfield(current_cycle, 'ankle_pos_FR2') && isfield(current_cycle, 'ankle_A_FR2') && isfield(current_cycle, 'ankle_b_FR2')
            transformed_pos_FR2 = zeros(size(current_cycle.ankle_pos_FR2));
            for t = 1:size(current_cycle.ankle_pos_FR2, 1)
                A_inv_FR2 = squeeze(current_cycle.ankle_A_FR2(t, :, :))';
                b_FR2 = current_cycle.ankle_b_FR2(t, :)';
                pos_FR2 = current_cycle.ankle_pos_FR2(t, :)';
                transformed_pos_FR2(t, :) = (A_inv_FR2 * (pos_FR2 - b_FR2))';
            end
            subplot(2, 3, 4);
            plot(current_cycle.time, transformed_pos_FR2, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), transformed_pos_FR2(1,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        if isfield(current_cycle, 'ankle_pos_FR2_velocity') && isfield(current_cycle, 'ankle_A_FR2')
            transformed_vel_FR2 = zeros(size(current_cycle.ankle_pos_FR2_velocity));
            for t = 1:size(current_cycle.ankle_pos_FR2_velocity, 1)
                A_inv_FR2 = squeeze(current_cycle.ankle_A_FR2(t, :, :))';
                vel_FR2 = current_cycle.ankle_pos_FR2_velocity(t, :)';
                transformed_vel_FR2(t, :) = (A_inv_FR2 * vel_FR2)';
            end
            subplot(2, 3, 5);
            plot(current_cycle.time, transformed_vel_FR2, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), transformed_vel_FR2(1,:), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        if isfield(current_cycle, 'ankle_orientation_FR2')
            subplot(2, 3, 6);
            plot(current_cycle.time, current_cycle.ankle_orientation_FR2, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_orientation_FR2(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
    end

    % Add origin marker and legends to all subplots for Figure 4
    for k = 1:6
        figure(4);
        subplot(2, 3, k);
        plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
%         legend('show', 'Location', 'best');
    end

    % Create figure for transformed phase space plots
    figure(5);
    set(gcf, 'Name', 'Transformed Phase Space Plots', 'Position', [300, 300, 1800, 800]);

    % Row 1: Original Frame (FR1) - Transformed
    subplot(2, 3, 1); hold on; title('Transformed Ankle Positions (FR1)'); xlabel('X Position'); ylabel('Y Position'); grid on; axis equal;
    subplot(2, 3, 2); hold on; title('Transformed Ankle Velocities (FR1)'); xlabel('X Velocity'); ylabel('Y Velocity'); grid on;
    subplot(2, 3, 3); hold on; title('Ankle Orientation (FR1)'); xlabel('Time (s)'); ylabel('Orientation (rad)'); grid on;

    % Row 2: Second Frame (FR2) - Transformed
    subplot(2, 3, 4); hold on; title('Transformed Ankle Positions (FR2)'); xlabel('X Position (FR2)'); ylabel('Y Position (FR2)'); grid on; axis equal;
    subplot(2, 3, 5); hold on; title('Transformed Ankle Velocities (FR2)'); xlabel('X Velocity (FR2)'); ylabel('Y Velocity (FR2)'); grid on;
    subplot(2, 3, 6); hold on; title('Ankle Orientation (FR2)'); xlabel('Time (s)'); ylabel('Orientation (rad)'); grid on;

    % Loop through all processed gait data to plot transformed data
    for i = 1:length(processed_gait_data)
        current_cycle = processed_gait_data{i};

        % --- Transformed Phase Space Plots ---
        figure(5);

        % FR1 Transformation and Plotting
        if isfield(current_cycle, 'ankle_pos_FR1') && isfield(current_cycle, 'ankle_A_FR1') && isfield(current_cycle, 'ankle_b_FR1')
            transformed_pos_FR1 = zeros(size(current_cycle.ankle_pos_FR1));
            for t = 1:size(current_cycle.ankle_pos_FR1, 1)
                A_inv_FR1 = squeeze(current_cycle.ankle_A_FR1(t, :, :))'; % Inverse of rotation matrix is its transpose
                b_FR1 = current_cycle.ankle_b_FR1(t, :)';
                pos_FR1 = current_cycle.ankle_pos_FR1(t, :)';
                transformed_pos_FR1(t, :) = (A_inv_FR1 * (pos_FR1 - b_FR1))';
            end
            subplot(2, 3, 1);
            plot(transformed_pos_FR1(:,1), transformed_pos_FR1(:,2), '--', 'HandleVisibility','off');
            plot(transformed_pos_FR1(1,1), transformed_pos_FR1(1,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        if isfield(current_cycle, 'ankle_pos_FR1_velocity') && isfield(current_cycle, 'ankle_A_FR1')
            transformed_vel_FR1 = zeros(size(current_cycle.ankle_pos_FR1_velocity));
            for t = 1:size(current_cycle.ankle_pos_FR1_velocity, 1)
                A_inv_FR1 = squeeze(current_cycle.ankle_A_FR1(t, :, :))';
                vel_FR1 = current_cycle.ankle_pos_FR1_velocity(t, :)';
                transformed_vel_FR1(t, :) = (A_inv_FR1 * vel_FR1)';
            end
            subplot(2, 3, 2);
            plot(transformed_vel_FR1(:,1), transformed_vel_FR1(:,2), '--', 'HandleVisibility','off');
            plot(transformed_vel_FR1(1,1), transformed_vel_FR1(1,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        if isfield(current_cycle, 'ankle_orientation_FR1')
            subplot(2, 3, 3);
            plot(current_cycle.time, current_cycle.ankle_orientation_FR1, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_orientation_FR1(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        % FR2 Transformation and Plotting
        if isfield(current_cycle, 'ankle_pos_FR2') && isfield(current_cycle, 'ankle_A_FR2') && isfield(current_cycle, 'ankle_b_FR2')
            transformed_pos_FR2 = zeros(size(current_cycle.ankle_pos_FR2));
            for t = 1:size(current_cycle.ankle_pos_FR2, 1)
                A_inv_FR2 = squeeze(current_cycle.ankle_A_FR2(t, :, :))';
                b_FR2 = current_cycle.ankle_b_FR2(t, :)';
                pos_FR2 = current_cycle.ankle_pos_FR2(t, :)';
                transformed_pos_FR2(t, :) = (A_inv_FR2 * (pos_FR2 - b_FR2))';
            end
            subplot(2, 3, 4);
            plot(transformed_pos_FR2(:,1), transformed_pos_FR2(:,2), '--', 'HandleVisibility','off');
            plot(transformed_pos_FR2(1,1), transformed_pos_FR2(1,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        if isfield(current_cycle, 'ankle_pos_FR2_velocity') && isfield(current_cycle, 'ankle_A_FR2')
            transformed_vel_FR2 = zeros(size(current_cycle.ankle_pos_FR2_velocity));
            for t = 1:size(current_cycle.ankle_pos_FR2_velocity, 1)
                A_inv_FR2 = squeeze(current_cycle.ankle_A_FR2(t, :, :))';
                vel_FR2 = current_cycle.ankle_pos_FR2_velocity(t, :)';
                transformed_vel_FR2(t, :) = (A_inv_FR2 * vel_FR2)';
            end
            subplot(2, 3, 5);
            plot(transformed_vel_FR2(:,1), transformed_vel_FR2(:,2), '--', 'HandleVisibility','off');
            plot(transformed_vel_FR2(1,1), transformed_vel_FR2(1,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end

        if isfield(current_cycle, 'ankle_orientation_FR2')
            subplot(2, 3, 6);
            plot(current_cycle.time, current_cycle.ankle_orientation_FR2, '--', 'HandleVisibility','off');
            plot(current_cycle.time(1), current_cycle.ankle_orientation_FR2(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Cycle %d Start', i));
        end
    end

    % Add origin marker and legends to all subplots for Figure 5
    for k = 1:6
        figure(5);
        subplot(2, 3, k);
        plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
%         legend('show', 'Location', 'best');
    end

    % Create figure for transformed phase space plots
    figure(5);
    set(gcf, 'Name', 'Transformed Phase Space Plots', 'Position', [300, 300, 1800, 800]);

    % Row 1: Original Frame (FR1) - Transformed
    subplot(2, 3, 1); hold on; title('Transformed Ankle Positions (FR1)'); xlabel('X Position'); ylabel('Y Position'); grid on; axis equal;
    subplot(2, 3, 2); hold on; title('Transformed Ankle Velocities (FR1)'); xlabel('X Velocity'); ylabel('Y Velocity'); grid on;
    subplot(2, 3, 3); hold on; title('Ankle Orientation (FR1)'); xlabel('Time (s)'); ylabel('Orientation (rad)'); grid on;

    % Row 2: Second Frame (FR2) - Transformed
    subplot(2, 3, 4); hold on; title('Transformed Ankle Positions (FR2)'); xlabel('X Position (FR2)'); ylabel('Y Position (FR2)'); grid on; axis equal;
    subplot(2, 3, 5); hold on; title('Transformed Ankle Velocities (FR2)'); xlabel('X Velocity (FR2)'); ylabel('Y Velocity (FR2)'); grid on;
    subplot(2, 3, 6); hold on; title('Ankle Orientation (FR2)'); xlabel('Time (s)'); ylabel('Orientation (rad)'); grid on;
end




