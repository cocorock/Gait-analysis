%% plot_gait_kinematics_v2: Plots raw/filtered positions, velocities, and accelerations for gait cycles.
%
% Credits:
%   Victor Ferman, Adrolab FEEC/UNICAMP
%   (Modified by Gemini)
%
% Description:
%   This function visualizes the processed gait data from apply_filtering_and_derivatives_v2.m.
%   It generates separate figures for right and left leg cycles, displaying raw and filtered
%   positions, velocities, and accelerations for each of the four main joints (right/left hip/knee).
%
% Input:
%   processed_data - struct: The processed gait data structure from apply_filtering_and_derivatives_v2.
%                     It must contain 'right_leg_cycles' and 'left_leg_cycles' fields,
%                     each being a struct array with raw, filtered, velocity, and acceleration data.
%   leg_to_plot    - (optional) string: 'right', 'left', or 'both'. Defaults to 'both'.
%
% Output:
%   Multiple figures displaying the kinematic data.

function plot_gait_kinematics_v2(processed_data, leg_to_plot)
    if nargin < 2
        leg_to_plot = 'both';
    end
    
    fprintf('\n=== PLOTTING GAIT KINEMATICS (V2) ===\n');
    
    time_standard = processed_data.time_standard;
    
    joint_names = {'Right Hip', 'Left Hip', 'Right Knee', 'Left Knee'};
    joint_fields_raw = {'right_hip_flex', 'left_hip_flex', 'right_knee_flex', 'left_knee_flex'};
    joint_fields_filtered = {'right_hip_flex_filtered', 'left_hip_flex_filtered', 'right_knee_flex_filtered', 'left_knee_flex_filtered'};
    joint_fields_velocity = {'right_hip_flex_velocity', 'left_hip_flex_velocity', 'right_knee_flex_velocity', 'left_knee_flex_velocity'};
    joint_fields_acceleration = {'right_hip_flex_acceleration', 'left_hip_flex_acceleration', 'right_knee_flex_acceleration', 'left_knee_flex_acceleration'};
    
    % --- Plot Right Leg Cycles ---
    if (strcmp(leg_to_plot, 'right') || strcmp(leg_to_plot, 'both')) && isfield(processed_data, 'right_leg_cycles') && ~isempty(processed_data.right_leg_cycles)
        fprintf('  Plotting Right Leg Cycles...\n');
        num_cycles = length(processed_data.right_leg_cycles);
        
        % Figure 1: Positions (Raw and Filtered)
        figure('Name', 'Right Leg Cycles - Positions', 'Position', [100, 100, 1200, 800]);
        sgtitle('Right Leg Cycles - Joint Positions (Raw vs Filtered)');
        for j = 1:length(joint_names)
            subplot(2, 2, j);
            hold on;
            for i = 1:num_cycles
                plot(time_standard, processed_data.right_leg_cycles(i).(joint_fields_raw{j}), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
                plot(time_standard, processed_data.right_leg_cycles(i).(joint_fields_filtered{j}), 'b-', 'LineWidth', 1.5);
            end
            hold off;
            title(joint_names{j});
            xlabel('Normalized Cycle Time');
            ylabel('Angle (degrees)');
            grid on;
            if j == 1
                legend('Raw', 'Filtered', 'Location', 'best');
            end
        end
        
        % Figure 2: Velocities
        figure('Name', 'Right Leg Cycles - Velocities', 'Position', [200, 200, 1200, 800]);
        sgtitle('Right Leg Cycles - Joint Velocities');
        for j = 1:length(joint_names)
            subplot(2, 2, j);
            hold on;
            for i = 1:num_cycles
                plot(time_standard, processed_data.right_leg_cycles(i).(joint_fields_velocity{j}), 'r-', 'LineWidth', 1.5);
            end
            hold off;
            title(joint_names{j});
            xlabel('Normalized Cycle Time');
            ylabel('Angular Velocity (deg/s)');
            grid on;
        end
        
        % Figure 3: Accelerations
        figure('Name', 'Right Leg Cycles - Accelerations', 'Position', [300, 300, 1200, 800]);
        sgtitle('Right Leg Cycles - Joint Accelerations');
        for j = 1:length(joint_names)
            subplot(2, 2, j);
            hold on;
            for i = 1:num_cycles
                plot(time_standard, processed_data.right_leg_cycles(i).(joint_fields_acceleration{j}), 'g-', 'LineWidth', 1.5);
            end
            hold off;
            title(joint_names{j});
            xlabel('Normalized Cycle Time');
            ylabel('Angular Acceleration (deg/s^2)');
            grid on;
        end
    else
        fprintf('  No right leg cycles to plot or plotting skipped.\n');
    end
    
    % --- Plot Left Leg Cycles ---
    if (strcmp(leg_to_plot, 'left') || strcmp(leg_to_plot, 'both')) && isfield(processed_data, 'left_leg_cycles') && ~isempty(processed_data.left_leg_cycles)
        fprintf('  Plotting Left Leg Cycles...\n');
        num_cycles = length(processed_data.left_leg_cycles);
        
        % Figure 4: Positions (Raw and Filtered)
        figure('Name', 'Left Leg Cycles - Positions', 'Position', [400, 400, 1200, 800]);
        sgtitle('Left Leg Cycles - Joint Positions (Raw vs Filtered)');
        for j = 1:length(joint_names)
            subplot(2, 2, j);
            hold on;
            for i = 1:num_cycles
                plot(time_standard, processed_data.left_leg_cycles(i).(joint_fields_raw{j}), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
                plot(time_standard, processed_data.left_leg_cycles(i).(joint_fields_filtered{j}), 'b-', 'LineWidth', 1.5);
            end
            hold off;
            title(joint_names{j});
            xlabel('Normalized Cycle Time');
            ylabel('Angle (degrees)');
            grid on;
            if j == 1
                legend('Raw', 'Filtered', 'Location', 'best');
            end
        end
        
        % Figure 5: Velocities
        figure('Name', 'Left Leg Cycles - Velocities', 'Position', [500, 500, 1200, 800]);
        sgtitle('Left Leg Cycles - Joint Velocities');
        for j = 1:length(joint_names)
            subplot(2, 2, j);
            hold on;
            for i = 1:num_cycles
                plot(time_standard, processed_data.left_leg_cycles(i).(joint_fields_velocity{j}), 'r-', 'LineWidth', 1.5);
            end
            hold off;
            title(joint_names{j});
            xlabel('Normalized Cycle Time');
            ylabel('Angular Velocity (deg/s)');
            grid on;
        end
        
        % Figure 6: Accelerations
        figure('Name', 'Left Leg Cycles - Accelerations', 'Position', [600, 600, 1200, 800]);
        sgtitle('Left Leg Cycles - Joint Accelerations');
        for j = 1:length(joint_names)
            subplot(2, 2, j);
            hold on;
            for i = 1:num_cycles
                plot(time_standard, processed_data.left_leg_cycles(i).(joint_fields_acceleration{j}), 'g-', 'LineWidth', 1.5);
            end
            hold off;
            title(joint_names{j});
            xlabel('Normalized Cycle Time');
            ylabel('Angular Acceleration (deg/s^2)');
            grid on;
        end
    else
        fprintf('  No left leg cycles to plot or plotting skipped.\n');
    end
    
    fprintf('Plotting complete!\n');
end
