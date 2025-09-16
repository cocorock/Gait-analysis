% create_filtering_plot.m

close all;
clc;
clear;

% Add necessary paths
addpath('../Functions_rev/');
addpath('../Gait Data/');

% --- Parameters ---
subject = '35';
filename = sprintf('../Gait Data/all_trajectories_ALL#%s.mat', subject);
sampling_freq = 120; % Hz
cutoff_freq = 6; % Hz for Butterworth filter
padding_duration = 0.2; % seconds of padding

% --- Load Data ---
load(filename, 'all_trajectories');

% Select a representative trajectory (e.g., the first one)
trajectory = all_trajectories{3};

% --- Extract a Single Gait Cycle ---
% Use left ankle position in FR1 to detect heel strikes
heel_strikes = detect_heel_strikes(trajectory.left_ankle_pos_FR1, sampling_freq);

if length(heel_strikes) < 2
    error('Not enough heel strikes detected to form a complete gait cycle.');
end

% Extract the first full gait cycle from the hip joint data
start_frame = heel_strikes(1);
end_frame = heel_strikes(2);
raw_hip_angle = trajectory.left_ankle_orientation_FR1(start_frame:end_frame);
time_vector = (0:length(raw_hip_angle)-1) / sampling_freq;

% --- Apply Circular Padding ---
padding_samples = round(padding_duration * sampling_freq);
padded_hip_angle = [raw_hip_angle(end-padding_samples+1:end); raw_hip_angle; raw_hip_angle(1:padding_samples)];

% Create a new time vector for the padded signal that starts at -padding_duration
padded_time_vector = (0:length(padded_hip_angle)-1) / sampling_freq - padding_duration;

% Adjust the original time vector to start at 0
time_vector = time_vector - time_vector(1);

% --- Apply Butterworth Filter ---
[b, a] = butter(4, cutoff_freq / (sampling_freq / 2), 'low');
filtered_padded_hip_angle = filtfilt(b, a, padded_hip_angle);

% Remove the padding
filtered_hip_angle = filtered_padded_hip_angle(padding_samples+1:end-padding_samples);

% --- Calculate Angular Velocity ---
angular_velocity = diff(filtered_hip_angle) * sampling_freq;
velocity_time_vector = time_vector(1:end-1);

% Calculate padded angular velocity
padded_angular_velocity = diff(filtered_padded_hip_angle) * sampling_freq;
padded_velocity_time_vector = padded_time_vector(1:end-1);

% --- Generate Plot ---
figure('Position', [100, 100, 800, 600]);

% Top Panel: Angular Position
ax1 = subplot(1,2, 1);
hold on;

% Plot the padded raw data
pre_pad_time_plot = -padding_duration;
post_pad_time_plot = time_vector(end) + padding_duration;

plot(padded_time_vector, padded_hip_angle, 'b--', 'LineWidth', 1, 'DisplayName', 'Padded Raw Data');

plot(time_vector, raw_hip_angle, 'b', 'LineWidth', 1.5, 'DisplayName', 'Raw Data');
plot(time_vector, filtered_hip_angle, 'r', 'LineWidth', 1.5, 'DisplayName', 'Filtered Data');

% Add vertical lines for cycle start and end
xline(0, 'k--', 'HandleVisibility', 'off', 'DisplayName', 'Gait Cycle Start/End');
xline(time_vector(end), 'k--', 'HandleVisibility', 'off');

% Shade the padded regions
x_patch_pre = [pre_pad_time_plot, 0, 0, pre_pad_time_plot];
x_patch_post = [time_vector(end), post_pad_time_plot, post_pad_time_plot, time_vector(end)];
ylim_vals = ylim;
patch(x_patch_pre, [ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Circular Padded Region');
patch(x_patch_post, [ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility','off');


title('Hip Joint Trajectory Processing');
xlabel('Time (s)');
ylabel('Angular Position (rad)');
legend('show');
grid on;

% Bottom Panel: Angular Velocity
ax2 = subplot(1, 2, 2);
hold on

plot(padded_velocity_time_vector, padded_angular_velocity, 'k--', 'LineWidth', 1, 'DisplayName', 'Padded Angular Velocity');
plot(velocity_time_vector, angular_velocity, 'k', 'LineWidth', 1.5, 'DisplayName', 'Angular Velocity');

xline(0, 'k--', 'HandleVisibility', 'off', 'DisplayName', 'Gait Cycle Start/End');
xline(velocity_time_vector(end), 'k--', 'HandleVisibility', 'off');

% Shade the padded regions for velocity plot
x_patch_pre_vel = [pre_pad_time_plot, 0, 0, pre_pad_time_plot];
x_patch_post_vel = [velocity_time_vector(end), post_pad_time_plot, post_pad_time_plot, velocity_time_vector(end)];
ylim_vals_vel = ylim;
patch(x_patch_pre_vel, [ylim_vals_vel(1) ylim_vals_vel(1) ylim_vals_vel(2) ylim_vals_vel(2)], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility','off');
patch(x_patch_post_vel, [ylim_vals_vel(1) ylim_vals_vel(1) ylim_vals_vel(2) ylim_vals_vel(2)], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility','off');

title('Angular Velocity from Filtered Data');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('show');
grid on;

% Align the x-axes
linkaxes([ax1, ax2], 'x');
xlim([pre_pad_time_plot, post_pad_time_plot]);

% --- Save the Plot ---
output_filename = '../Filter/trajectory_filtering_plot.png';
saveas(gcf, output_filename);
fprintf('Plot saved to %s\n', output_filename);