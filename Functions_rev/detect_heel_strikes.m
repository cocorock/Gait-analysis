% Functions_rev/detect_heel_strikes.m
function heel_strike_indices = detect_heel_strikes(ankle_pos_FR2_trajectory, sample_rate)
    % DETECT_HEEL_STRIKES Detects robust heel strike events from ankle position trajectory
    %   in a specific frame of reference (FR2). Heel strike is identified as
    %   the moment of largest change (peak in acceleration magnitude) in the trajectory.
    %
    %   Inputs:
    %     ankle_pos_FR2_trajectory: Nx2 matrix where N is the number of frames,
    %                               column 1 is X-coordinate, column 2 is Y-coordinate
    %                               in the FR2 frame of reference.
    %     sample_rate: Sampling rate of the data in Hz.
    %
    %   Output:
    %     heel_strike_indices: Column vector of indices where heel strikes occur.

    if nargin < 2
        error('Not enough input arguments. Usage: detect_heel_strikes(ankle_pos_FR2_trajectory, sample_rate)');
    end

    % --- Debugging Flag ---
    show_debug_plot = true; % Set to true to show plots for debugging

    % --- Parameters for robust detection ---
    window_size = 5; % Smoothing window size (frames)
    min_gait_duration = 0.75; % Minimum expected gait cycle duration (seconds)
    min_frames_between_cycles = round(min_gait_duration * sample_rate);
    
    % Time step for velocity/acceleration calculation
    dt = 1 / sample_rate;

    % 1. Calculate velocity of the trajectory
    velocity_trajectory = calculate_velocity(ankle_pos_FR2_trajectory, dt);

    % --- IMPORTANT: Remove the last element due to circular reference issue ---
    % The last element of velocity_trajectory can be erroneous when the input
    % trajectory is not a perfect cycle.
    velocity_trajectory = velocity_trajectory(1:end-1, :);

    % 2. Calculate acceleration from velocity
    acceleration_trajectory = calculate_velocity(velocity_trajectory, dt); % Re-using calculate_velocity for derivative

    % --- IMPORTANT: Remove the last element due to circular reference issue (again for acceleration) ---
    acceleration_trajectory = acceleration_trajectory(1:end-1, :);

    % 3. Calculate the magnitude of the acceleration
    acceleration_magnitude = sqrt(sum(acceleration_trajectory.^2, 2));

    % 4. Smooth the acceleration magnitude
    acceleration_smooth = smooth(acceleration_magnitude, window_size);

    % 5. Detect peaks in the smoothed acceleration magnitude
    % Heel strike is often associated with a peak in acceleration.
    
    % Calculate prominence threshold based on data range
    accel_range = max(acceleration_smooth) - min(acceleration_smooth);
    min_prominence_threshold = accel_range * 0.10; % 10% of the acceleration range as minimum prominence

    [pks, heel_strike_indices] = findpeaks(acceleration_smooth, ...
        'MinPeakDistance', min_frames_between_cycles, ...
        'MinPeakProminence', min_prominence_threshold);

    % Ensure indices are column vector
    heel_strike_indices = heel_strike_indices(:);

    % --- Debugging Plots ---
    if show_debug_plot
        figure('Name', 'Heel Strike Detection Debug', 'Position', [100, 100, 1000, 700]);

        % Plot 1: Original Ankle Y-trajectory (for context, though not directly used for detection anymore)
        subplot(3,1,1);
        plot(ankle_pos_FR2_trajectory(:,2), 'b');
        title('Original Ankle Y-trajectory (FR2)');
        xlabel('Frame');
        ylabel('Y-position');
        grid on;

        % Plot 2: Acceleration Magnitude
        subplot(3,1,2);
        plot(acceleration_magnitude, 'k');
        hold on;
        plot(acceleration_smooth, 'r', 'LineWidth', 1.5);
        title('Acceleration Magnitude (Black: Original, Red: Smoothed)');
        xlabel('Frame');
        ylabel('Acceleration');
        grid on;
        legend('Original Acceleration', 'Smoothed Acceleration', 'Location', 'best');

        % Plot 3: Smoothed Acceleration Magnitude with Detected Heel Strikes
        subplot(3,1,3);
        plot(acceleration_smooth, 'r', 'LineWidth', 1.5);
        hold on;
        if ~isempty(heel_strike_indices)
            plot(heel_strike_indices, acceleration_smooth(heel_strike_indices), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
            text(heel_strike_indices, acceleration_smooth(heel_strike_indices), num2str(heel_strike_indices), ...
                 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'g');
        end
        title(sprintf('Smoothed Acceleration Magnitude with Detected Heel Strikes (%d found)', length(heel_strike_indices)));
        xlabel('Frame');
        ylabel('Acceleration');
        grid on;
        if ~isempty(heel_strike_indices)
            legend('Smoothed Acceleration', 'Detected Heel Strikes', 'Location', 'best');
        else
            legend('Smoothed Acceleration', 'Location', 'best');
        end
        
        % Add a line for the prominence threshold
        yline(min(acceleration_smooth) + min_prominence_threshold, ':', 'Prominence Threshold');
    end
end