function plot_root_motion(D)
% Plots root position and rotation from AMC motion data
% Input: D - matrix from amc_to_matrix function (N x 62)
% No rotation applied, only scaling

% Extract and scale root position (columns 1-3)
root_pos = D(:, 1:3) * 2.223 * 25.4 / 1000;  % Position in meters
root_rot = D(:, 4:6);  % Rotation data (degrees)

% Create time vector (frame numbers)
frames = 1:size(D, 1);

% Create figure with subplots
figure('Position', [100, 100, 1200, 800]);

% 2D Plot 1: Position components over time
subplot(2, 3, 1);
plot(frames, root_pos(:,1), 'r-', 'LineWidth', 1.5); hold on;
plot(frames, root_pos(:,2), 'g-', 'LineWidth', 1.5);
plot(frames, root_pos(:,3), 'b-', 'LineWidth', 1.5);
xlabel('Frame');
ylabel('Position (m)');
title('Root Position Over Time');
legend('X', 'Y', 'Z', 'Location', 'best');
grid on;

% 2D Plot 2: Rotation components over time
subplot(2, 3, 2);
plot(frames, root_rot(:,1), 'r-', 'LineWidth', 1.5); hold on;
plot(frames, root_rot(:,2), 'g-', 'LineWidth', 1.5);
plot(frames, root_rot(:,3), 'b-', 'LineWidth', 1.5);
xlabel('Frame');
ylabel('Rotation (degrees)');
title('Root Rotation Over Time');
legend('RX', 'RY', 'RZ', 'Location', 'best');
grid on;

% 2D Plot 3: X-Y trajectory (top view)
subplot(2, 3, 3);
plot(root_pos(:,1), root_pos(:,2), 'b-', 'LineWidth', 1.5);
hold on;
plot(root_pos(1,1), root_pos(1,2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start
plot(root_pos(end,1), root_pos(end,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % End
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Root Trajectory (Top View)');
legend('Path', 'Start', 'End', 'Location', 'best');
grid on;

% 2D Plot 4: Z-X trajectory (side view, Z horizontal)
subplot(2, 3, 4);
plot(root_pos(:,3), root_pos(:,1), 'b-', 'LineWidth', 1.5);
hold on;
plot(root_pos(1,3), root_pos(1,1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start
plot(root_pos(end,3), root_pos(end,1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % End
xlabel('Z Position (m)');
ylabel('X Position (m)');
title('Root Trajectory (Side View, Z Horizontal)');
legend('Path', 'Start', 'End', 'Location', 'best');
grid on;

% 2D Plot 5: Z-Y trajectory (front view, Z horizontal)
subplot(2, 3, 5);
plot(root_pos(:,3), root_pos(:,2), 'b-', 'LineWidth', 1.5);
hold on;
plot(root_pos(1,3), root_pos(1,2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start
plot(root_pos(end,3), root_pos(end,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % End
xlabel('Z Position (m)');
ylabel('Y Position (m)');
title('Root Trajectory (Front View, Z Horizontal)');
legend('Path', 'Start', 'End', 'Location', 'best');
grid on;

% 3D Plot: Complete 3D trajectory, Z horizontal
subplot(2, 3, 6);
plot3(root_pos(:,1), root_pos(:,2), root_pos(:,3), 'b-', 'LineWidth', 2);
hold on;
plot3(root_pos(1,1), root_pos(1,2), root_pos(1,3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); % Start
plot3(root_pos(end,1), root_pos(end,2), root_pos(end,3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % End
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('Root 3D Trajectory (Z Horizontal)');
legend('Path', 'Start', 'End', 'Location', 'best');
grid on;
% axis equal;
view([1 0 0]); % Z axis horizontal (points to the right)

% Add overall title
sgtitle('Root Motion Analysis (Scaled, No Rotation)', 'FontSize', 16, 'FontWeight', 'bold');


figure()
plot3(root_pos(:,1), root_pos(:,2), root_pos(:,3), 'b-', 'LineWidth', 2);
hold on;
plot3(root_pos(1,1), root_pos(1,2), root_pos(1,3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); % Start
plot3(root_pos(end,1), root_pos(end,2), root_pos(end,3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % End
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('Root 3D Trajectory (Z Horizontal)');
legend('Path', 'Start', 'End', 'Location', 'best');
grid on;
axis equal;
view([0 0 -1]); % Z axis horizontal (points to the right)
ylim([0 1.2])
xlim([-0.5 0.5])
end