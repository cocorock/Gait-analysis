% Functions_rev/detect_heel_strikes.m
function heel_strike_indices = detect_heel_strikes(ankle_pos_trajectory, sample_rate)
    % DETECT_HEEL_STRIKES Detects robust heel strike events from ankle position trajectory.
    %   heel_strike_indices = DETECT_HEEL_STRIKES(ankle_pos_trajectory, sample_rate)
    %   finds the indices corresponding to heel strike events, which are
    %   assumed to be the lowest points (local minima) in the vertical (Y)
    %   component of the ankle trajectory. This version incorporates smoothing
    %   and uses findpeaks for more robust minima detection.
    %
    %   Inputs:
    %     ankle_pos_trajectory: Nx2 matrix where N is the number of frames,
    %                           column 1 is X-coordinate, column 2 is Y-coordinate.
    %     sample_rate: Sampling rate of the data in Hz.
    %
    %   Output:
    %     heel_strike_indices: Column vector of indices where heel strikes occur.

    if nargin < 2
        error('Not enough input arguments. Usage: detect_heel_strikes(ankle_pos_trajectory, sample_rate)');
    end

    y_trajectory = ankle_pos_trajectory(:, 2); % Get the Y-coordinate

    % --- Parameters for robust detection ---
    window_size = 5; % Smoothing window size (frames)
    min_gait_duration = 0.75; % Minimum expected gait cycle duration (seconds)
    min_frames_between_cycles = round(min_gait_duration * sample_rate);

    % Calculate prominence threshold based on data range
    y_range = max(y_trajectory) - min(y_trajectory);
    min_prominence_threshold = y_range * 0.10; % 10% of the Y-range as minimum prominence

    % 1. Smooth the signal
    y_smooth = smooth(y_trajectory, window_size);

    % 2. Detect minima using findpeaks on the inverted signal
    % findpeaks finds peaks, so we invert the signal to find minima
    [~, heel_strike_indices] = findpeaks(-y_smooth, ...
        'MinPeakDistance', min_frames_between_cycles, ...
        'MinPeakProminence', min_prominence_threshold);

    % Ensure indices are column vector
    heel_strike_indices = heel_strike_indices(:);
end