function [gait_cycles, heel_strikes] = segment_gait_cycles(filename)
% SEGMENT_GAIT_CYCLES - Detects heel strikes and segments gait cycles from AMC data
%
% Input:
%   filename - path to the AMC file
%
% Output:
%   gait_cycles - struct array containing segmented gait cycles
%   heel_strikes - struct with detected heel strike frames
%
% Sample rate: 120 Hz

    % Read and parse AMC file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    % Skip header lines
    line = fgetl(fid);
    while ischar(line) && ~strcmp(strtrim(line), ':DEGREES')
        line = fgetl(fid);
    end

    % Initialize data storage
    frame_data = [];
    current_frame = 0;

    % Read frame data
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(strtrim(line))
            tokens = strsplit(strtrim(line));

            % Check if this is a frame number line
            if length(tokens) == 1 && ~isnan(str2double(tokens{1}))
                current_frame = str2double(tokens{1});
            else
                % This is joint data
                if current_frame > 0 && length(tokens) >= 2
                    joint_name = tokens{1};
                    values = str2double(tokens(2:end));

                    % Store heel marker data (assuming heel markers exist)
                    % We'll use foot segment orientations as proxy for heel position
                    if strcmp(joint_name, 'rfoot') || strcmp(joint_name, 'lfoot')
                        frame_data(current_frame).(joint_name) = values;
                    end

                    % Also store root position for reference
                    if strcmp(joint_name, 'root')
                        frame_data(current_frame).root = values;
                    end
                end
            end
        end
    end
    fclose(fid);

    % Extract time vector (120 Hz sampling rate)
    num_frames = length(frame_data);
    time = (0:num_frames-1) / 120; % 120 Hz sample rate

    % Extract heel trajectories
    % Since AMC doesn't have direct heel markers, we'll use foot segment data
    % and estimate heel position from foot orientation and root position

    right_heel_y = zeros(num_frames, 1);
    left_heel_y = zeros(num_frames, 1);
    right_heel_vel = zeros(num_frames, 1);
    left_heel_vel = zeros(num_frames, 1);

    for i = 1:num_frames
        if isfield(frame_data(i), 'root') && isfield(frame_data(i), 'rfoot') && isfield(frame_data(i), 'lfoot')
            % Estimate heel height from root position and foot orientation
            root_y = frame_data(i).root(2); % Y position of root

            % Use foot pitch angle to estimate heel height relative to root
            if length(frame_data(i).rfoot) >= 2
                rfoot_pitch = frame_data(i).rfoot(2); % Foot pitch angle
                right_heel_y(i) = root_y - 100 * sind(rfoot_pitch); % Approximate heel height
            end

            if length(frame_data(i).lfoot) >= 2
                lfoot_pitch = frame_data(i).lfoot(2); % Foot pitch angle
                left_heel_y(i) = root_y - 100 * sind(lfoot_pitch); % Approximate heel height
            end
        end
    end

    % Calculate heel velocities (forward velocity)
    if num_frames > 1
        dt = 1/120; % Time step
        right_heel_vel(2:end) = diff(right_heel_y) / dt;
        left_heel_vel(2:end) = diff(left_heel_y) / dt;

        % Smooth velocities
        if num_frames > 10
            right_heel_vel = smooth(right_heel_vel, 5);
            left_heel_vel = smooth(left_heel_vel, 5);
        end
    end

    % Detect heel strikes using local minima in heel height
    % Combined with low forward velocity

    % Parameters for heel strike detection
    min_cycle_duration = 0.8; % Minimum gait cycle duration (seconds)
    min_frames_between_strikes = round(min_cycle_duration * 120);

    % Detect right heel strikes
    [~, right_minima] = findpeaks(-right_heel_y, 'MinPeakDistance', min_frames_between_strikes, ...
                                  'MinPeakHeight', -max(right_heel_y) * 0.1);

    % Detect left heel strikes
    [~, left_minima] = findpeaks(-left_heel_y, 'MinPeakDistance', min_frames_between_strikes, ...
                                 'MinPeakHeight', -max(left_heel_y) * 0.1);

    % Refine heel strikes by checking velocity
    right_heel_strikes = [];
    for i = 1:length(right_minima)
        idx = right_minima(i);
        if idx > 5 && idx < num_frames - 5
            % Check if velocity is low around this point
            vel_window = abs(right_heel_vel(idx-2:idx+2));
            if mean(vel_window) < std(right_heel_vel) * 0.5
                right_heel_strikes = [right_heel_strikes, idx];
            end
        end
    end

    left_heel_strikes = [];
    for i = 1:length(left_minima)
        idx = left_minima(i);
        if idx > 5 && idx < num_frames - 5
            % Check if velocity is low around this point
            vel_window = abs(left_heel_vel(idx-2:idx+2));
            if mean(vel_window) < std(left_heel_vel) * 0.5
                left_heel_strikes = [left_heel_strikes, idx];
            end
        end
    end

    % Store heel strike information
    heel_strikes.right_frames = right_heel_strikes;
    heel_strikes.left_frames = left_heel_strikes;
    heel_strikes.right_times = right_heel_strikes / 120;
    heel_strikes.left_times = left_heel_strikes / 120;

    % Segment gait cycles
    gait_cycles = [];
    cycle_count = 0;

    % Right leg gait cycles
    for i = 1:length(right_heel_strikes)-1
        cycle_count = cycle_count + 1;
        start_frame = right_heel_strikes(i);
        end_frame = right_heel_strikes(i+1);

        gait_cycles(cycle_count).leg = 'right';
        gait_cycles(cycle_count).start_frame = start_frame;
        gait_cycles(cycle_count).end_frame = end_frame;
        gait_cycles(cycle_count).start_time = start_frame / 120;
        gait_cycles(cycle_count).end_time = end_frame / 120;
        gait_cycles(cycle_count).duration = (end_frame - start_frame) / 120;
        gait_cycles(cycle_count).frames = start_frame:end_frame;
        gait_cycles(cycle_count).time_normalized = linspace(0, 100, end_frame - start_frame + 1);
    end

    % Left leg gait cycles
    for i = 1:length(left_heel_strikes)-1
        cycle_count = cycle_count + 1;
        start_frame = left_heel_strikes(i);
        end_frame = left_heel_strikes(i+1);

        gait_cycles(cycle_count).leg = 'left';
        gait_cycles(cycle_count).start_frame = start_frame;
        gait_cycles(cycle_count).end_frame = end_frame;
        gait_cycles(cycle_count).start_time = start_frame / 120;
        gait_cycles(cycle_count).end_time = end_frame / 120;
        gait_cycles(cycle_count).duration = (end_frame - start_frame) / 120;
        gait_cycles(cycle_count).frames = start_frame:end_frame;
        gait_cycles(cycle_count).time_normalized = linspace(0, 100, end_frame - start_frame + 1);
    end

    % Plot heel strike detection results
    figure('Name', 'Heel Strike Detection', 'Position', [100, 100, 1200, 800]);

    subplot(3,1,1);
    plot(time, right_heel_y, 'r-', 'LineWidth', 1.5);
    hold on;
    plot(time, left_heel_y, 'b-', 'LineWidth', 1.5);
    if ~isempty(right_heel_strikes)
        plot(time(right_heel_strikes), right_heel_y(right_heel_strikes), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
    if ~isempty(left_heel_strikes)
        plot(time(left_heel_strikes), left_heel_y(left_heel_strikes), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    end
    xlabel('Time (s)');
    ylabel('Heel Height (mm)');
    title('Heel Height and Detected Heel Strikes');
    legend('Right Heel', 'Left Heel', 'Right Heel Strikes', 'Left Heel Strikes');
    grid on;

    subplot(3,1,2);
    plot(time, right_heel_vel, 'r-', 'LineWidth', 1);
    hold on;
    plot(time, left_heel_vel, 'b-', 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Heel Velocity (mm/s)');
    title('Heel Vertical Velocity');
    legend('Right Heel', 'Left Heel');
    grid on;

    subplot(3,1,3);
    % Plot gait cycle durations
    if ~isempty(gait_cycles)
        right_cycles = find(strcmp({gait_cycles.leg}, 'right'));
        left_cycles = find(strcmp({gait_cycles.leg}, 'left'));

        if ~isempty(right_cycles)
            bar(right_cycles, [gait_cycles(right_cycles).duration], 'r', 'FaceAlpha', 0.7);
            hold on;
        end
        if ~isempty(left_cycles)
            bar(left_cycles, [gait_cycles(left_cycles).duration], 'b', 'FaceAlpha', 0.7);
        end

        xlabel('Gait Cycle Number');
        ylabel('Duration (s)');
        title('Gait Cycle Durations');
        legend('Right Leg', 'Left Leg');
        grid on;
    end

    % Display summary
    fprintf('\n=== GAIT CYCLE SEGMENTATION RESULTS ===\n');
    fprintf('Total frames analyzed: %d\n', num_frames);
    fprintf('Total duration: %.2f seconds\n', num_frames/120);
    fprintf('Right heel strikes detected: %d\n', length(right_heel_strikes));
    fprintf('Left heel strikes detected: %d\n', length(left_heel_strikes));
    fprintf('Total gait cycles identified: %d\n', length(gait_cycles));

    if ~isempty(gait_cycles)
        right_cycles = find(strcmp({gait_cycles.leg}, 'right'));
        left_cycles = find(strcmp({gait_cycles.leg}, 'left'));

        if ~isempty(right_cycles)
            fprintf('Right leg cycles: %d (avg duration: %.2f ± %.2f s)\n', ...
                    length(right_cycles), ...
                    mean([gait_cycles(right_cycles).duration]), ...
                    std([gait_cycles(right_cycles).duration]));
        end

        if ~isempty(left_cycles)
            fprintf('Left leg cycles: %d (avg duration: %.2f ± %.2f s)\n', ...
                    length(left_cycles), ...
                    mean([gait_cycles(left_cycles).duration]), ...
                    std([gait_cycles(left_cycles).duration]));
        end
    end

end