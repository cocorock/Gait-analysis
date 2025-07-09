
% Functions_rev/detect_heel_strikes.m
function heel_strike_indices = detect_heel_strikes(ankle_pos_FR2_trajectory, sample_rate)
    % DETECT_HEEL_STRIKES Detects robust heel strike events from ankle position trajectory
    %   in a specific frame of reference (FR2). Heel strike is identified as
    %   the moment of largest change (peak in speed magnitude) in the trajectory.
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
    
    % Time step for velocity calculation
    dt = 1 / sample_rate;

    % 1. Calculate velocity of the trajectory
    % Assuming calculate_velocity function is available in the path
    velocity_trajectory = calculate_velocity(ankle_pos_FR2_trajectory, dt);

    % --- IMPORTANT: Remove the last element due to circular reference issue ---
    % The last element of velocity_trajectory can be erroneous when the input
    % trajectory is not a perfect cycle.
    velocity_trajectory = velocity_trajectory(1:end-1, :);

    % 2. Calculate the magnitude of the velocity (speed)
    speed_magnitude = sqrt(sum(velocity_trajectory.^2, 2));

    % 3. Smooth the speed magnitude
    speed_smooth = smooth(speed_magnitude, window_size);

    % 4. Detect peaks in the smoothed speed magnitude
    % Heel strike is often associated with a peak in speed or a sudden change.
    % We are looking for the "largest change", which implies a peak in speed.
    
    % Calculate prominence threshold based on data range
    speed_range = max(speed_smooth) - min(speed_smooth);
    min_prominence_threshold = speed_range * 0.10; % 10% of the speed range as minimum prominence

    [pks, heel_strike_indices] = findpeaks(speed_smooth, ...
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

        % Plot 2: Speed Magnitude
        subplot(3,1,2);
        plot(speed_magnitude, 'k');
        hold on;
        plot(speed_smooth, 'r', 'LineWidth', 1.5);
        title('Speed Magnitude (Black: Original, Red: Smoothed)');
        xlabel('Frame');
        ylabel('Speed');
        grid on;
        legend('Original Speed', 'Smoothed Speed', 'Location', 'best');

        % Plot 3: Smoothed Speed Magnitude with Detected Heel Strikes
        subplot(3,1,3);
        plot(speed_smooth, 'r', 'LineWidth', 1.5);
        hold on;
        if ~isempty(heel_strike_indices)
            plot(heel_strike_indices, speed_smooth(heel_strike_indices), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
            text(heel_strike_indices, speed_smooth(heel_strike_indices), num2str(heel_strike_indices), ...
                 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'g');
        end
        title(sprintf('Smoothed Speed Magnitude with Detected Heel Strikes (%d found)', length(heel_strike_indices)));
        xlabel('Frame');
        ylabel('Speed');
        grid on;
        if ~isempty(heel_strike_indices)
            legend('Smoothed Speed', 'Detected Heel Strikes', 'Location', 'best');
        else
            legend('Smoothed Speed', 'Location', 'best');
        end
        
        % Add a line for the prominence threshold
        yline(min(speed_smooth) + min_prominence_threshold, ':', 'Prominence Threshold');
    end
end
