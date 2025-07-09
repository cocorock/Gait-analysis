function trajectories = plot_root_and_ankles_trajectory(asf_file, amc_file, plot_figures)
% Plots the trajectory of the root, knees, and ankles.
% Calculates and plots the ankle orientation.
% Inputs:
%   asf_file - path to the .asf file
%   amc_file - path to the .amc file
%   plot_figures - flag to plot figures or not
% Outputs:
%   trajectories - struct containing trajectory and orientation data

    % Read bone lengths from ASF file
    bone_lengths = read_asf_lengths(asf_file);

    % Read motion data from AMC file
    D = amc_to_matrix(amc_file);
    frames = 0:1/120:size(D,1)/120; 
    frames = frames(1:end-1)';
    
    % Extract root position and rotation
    root_pos = D(:, 1:3) * 1/0.45 * 25.4 / 1000; % Scale to meters
    root_rot = D(:, 4:6);

    % Get joint angles
    lfemur_rot = -D(:, 56); %Flexion/Extesion  %56:58
    rfemur_rot = -D(:, 49); %Flexion/Extesion  %49:51
    ltibia_rot = -D(:, 59);   %Flexion/Extesion 
    rtibia_rot = -D(:, 52);   %Flexion/Extesion 

    % Initialize position matrices
    left_knee_pos = zeros(size(root_pos));
    right_knee_pos = zeros(size(root_pos));
    left_ankle_pos = zeros(size(root_pos));
    right_ankle_pos = zeros(size(root_pos));
    left_hip_pos_all = zeros(size(root_pos));
    right_hip_pos_all = zeros(size(root_pos));

    for i = 1:size(D, 1)
        % Get root position and rotation for the current frame
        hip_pos = root_pos(i, :) * [1 0 0; 0 1 0; 0 0 1];
        rx = root_rot(i, 1);
        ry = 0;%root_rot(i, 2);
        rz = 0;%root_rot(i, 3);

        % Create rotation matrices
        R_x = [1, 0, 0; 0, cosd(rx), -sind(rx); 0, sind(rx), cosd(rx)];
        R_y = [cosd(ry), 0, sind(ry); 0, 1, 0; -sind(ry), 0, cosd(ry)];
        R_z = [cosd(rz), -sind(rz), 0; sind(rz), cosd(rz), 0; 0, 0, 1];
        R = R_z * R_y * R_x; % Combined rotation matrix

        % Left leg kinematics
        left_hip_offset = [-bone_lengths.lhipjoint, 0, 0]';
        rotated_left_hip_offset = R * left_hip_offset;
        left_hip_pos = hip_pos + rotated_left_hip_offset';
        left_hip_pos_all(i, :) = left_hip_pos;

        knee_vec_l = bone_lengths.lfemur * [0, -cosd(lfemur_rot(i)), sind(lfemur_rot(i))]';
        rotated_knee_vec_l = R * knee_vec_l;
        left_knee_pos(i, :) = left_hip_pos + rotated_knee_vec_l';

        ankle_vec_l = bone_lengths.ltibia * [0, -cosd(lfemur_rot(i) + ltibia_rot(i)), sind(lfemur_rot(i) + ltibia_rot(i))]';
        rotated_ankle_vec_l = R * ankle_vec_l;
        left_ankle_pos(i, :) = left_knee_pos(i, :) + rotated_ankle_vec_l';

        % Right leg kinematics
        right_hip_offset = [bone_lengths.rhipjoint, 0, 0]';
        rotated_right_hip_offset = R * right_hip_offset;
        right_hip_pos = hip_pos + rotated_right_hip_offset';
        right_hip_pos_all(i, :) = right_hip_pos;

        knee_vec_r = bone_lengths.rfemur * [0, -cosd(rfemur_rot(i)), sind(rfemur_rot(i))]';
        rotated_knee_vec_r = R * knee_vec_r;
        right_knee_pos(i, :) = right_hip_pos + rotated_knee_vec_r';

        ankle_vec_r = bone_lengths.rtibia * [0, -cosd(rfemur_rot(i) + rtibia_rot(i)), sind(rfemur_rot(i) + rtibia_rot(i))]';
        rotated_ankle_vec_r = R * ankle_vec_r;
        right_ankle_pos(i, :) = right_knee_pos(i, :) + rotated_ankle_vec_r';
    end

    %% second frame of reference (hip)
    left_knee_pos_FR2 = zeros(size(root_pos));
    right_knee_pos_FR2 = zeros(size(root_pos));
    left_ankle_pos_FR2 = zeros(size(root_pos));
    right_ankle_pos_FR2 = zeros(size(root_pos));
    left_hip_pos_all_FR2 = zeros(size(root_pos));
    right_hip_pos_all_FR2 = zeros(size(root_pos));

    for i = 1:size(D, 1)
        hip_pos = 0; 
        rx = 0; 
        ry = 0;
        rz = 0;

        R_x = [1, 0, 0; 0, cosd(rx), -sind(rx); 0, sind(rx), cosd(rx)];
        R_y = [cosd(ry), 0, sind(ry); 0, 1, 0; -sind(ry), 0, cosd(ry)];
        R_z = [cosd(rz), -sind(rz), 0; sind(rz), cosd(rz), 0; 0, 0, 1];
        R = R_z * R_y * R_x; 

        left_hip_offset = [-bone_lengths.lhipjoint, 0, 0]';
        rotated_left_hip_offset = R * left_hip_offset;
        left_hip_pos = hip_pos + rotated_left_hip_offset';
        left_hip_pos_all_FR2(i, :) = left_hip_pos;

        knee_vec_l = bone_lengths.lfemur * [0, -cosd(lfemur_rot(i)), sind(lfemur_rot(i))]';
        rotated_knee_vec_l = R * knee_vec_l;
        left_knee_pos_FR2(i, :) = left_hip_pos + rotated_knee_vec_l';

        ankle_vec_l = bone_lengths.ltibia * [0, -cosd(lfemur_rot(i) + ltibia_rot(i)), sind(lfemur_rot(i) + ltibia_rot(i))]';
        rotated_ankle_vec_l = R * ankle_vec_l;
        left_ankle_pos_FR2(i, :) = left_knee_pos_FR2(i, :) + rotated_ankle_vec_l';

        right_hip_offset = [bone_lengths.rhipjoint, 0, 0]';
        rotated_right_hip_offset = R * right_hip_offset;
        right_hip_pos = hip_pos + rotated_right_hip_offset';
        right_hip_pos_all_FR2(i, :) = right_hip_pos;

        knee_vec_r = bone_lengths.rfemur * [0, -cosd(rfemur_rot(i)), sind(rfemur_rot(i))]';
        rotated_knee_vec_r = R * knee_vec_r;
        right_knee_pos_FR2(i, :) = right_hip_pos + rotated_knee_vec_r';

        ankle_vec_r = bone_lengths.rtibia * [0, -cosd(rfemur_rot(i) + rtibia_rot(i)), sind(rfemur_rot(i) + rtibia_rot(i))]';
        rotated_ankle_vec_r = R * ankle_vec_r;
        right_ankle_pos_FR2(i, :) = right_knee_pos_FR2(i, :) + rotated_ankle_vec_r';
    end

    % Calculate ankle orientation
    tibia_vec_l = left_ankle_pos - left_knee_pos;
    tibia_vec_r = right_ankle_pos - right_knee_pos;
    left_tibia_angle = atan2(tibia_vec_l(:,2), tibia_vec_l(:,3));
    right_tibia_angle = atan2(tibia_vec_r(:,2), tibia_vec_r(:,3));
    left_ankle_orientation = left_tibia_angle + pi/2;
    right_ankle_orientation = right_tibia_angle + pi/2;

    tibia_vec_l_FR2 = left_ankle_pos_FR2 - left_knee_pos_FR2;
    tibia_vec_r_FR2 = right_ankle_pos_FR2 - right_knee_pos_FR2;
    left_tibia_angle_FR2 = atan2(tibia_vec_l_FR2(:,2), tibia_vec_l_FR2(:,3));
    right_tibia_angle_FR2 = atan2(tibia_vec_r_FR2(:,2), tibia_vec_r_FR2(:,3));
    left_ankle_orientation_FR2 = left_tibia_angle_FR2 + pi/2;
    right_ankle_orientation_FR2 = right_tibia_angle_FR2 + pi/2;

    if plot_figures
        
        inter_space = 12;

        figure('Position', [100, 100, 1200, 800]);

        subplot(2, 3, 1);
        plot(frames, root_pos(:,1), 'k-', 'LineWidth', 1.5); hold on;
        plot(frames, root_pos(:,2), 'k-', 'LineWidth', 1.5);
        plot(frames, root_pos(:,3), 'k-', 'LineWidth', 1.5);
        plot(frames, left_ankle_pos(:,1), 'r-', 'LineWidth', 1);
        plot(frames, left_ankle_pos(:,2), 'r-', 'LineWidth', 1);
        plot(frames, left_ankle_pos(:,3), 'r-', 'LineWidth', 1);
        plot(frames, right_ankle_pos(:,1), 'b-', 'LineWidth', 1);
        plot(frames, right_ankle_pos(:,2), 'b-', 'LineWidth', 1);
        plot(frames, right_ankle_pos(:,3), 'b-', 'LineWidth', 1);
        plot(frames, left_knee_pos(:,1), 'r:', 'LineWidth', 1.5);
        plot(frames, left_knee_pos(:,2), 'r:', 'LineWidth', 1.5);
        plot(frames, left_knee_pos(:,3), 'r:', 'LineWidth', 1.5);
        plot(frames, right_knee_pos(:,1), 'b:', 'LineWidth', 1.5);
        plot(frames, right_knee_pos(:,2), 'b:', 'LineWidth', 1.5);
        plot(frames, right_knee_pos(:,3), 'b:', 'LineWidth', 1.5);
        xlabel('Time');
        ylabel('Position (m)');
        title('Position Over Time');
        legend('Root X', 'Root Y', 'Root Z', 'L Ankle X', 'L Ankle Y', 'L Ankle Z', 'R Ankle X', 'R Ankle Y', 'R Ankle Z', 'L Knee X', 'L Knee Y', 'L Knee Z', 'R Knee X', 'R Knee Y', 'R Knee Z', 'Location', 'best');
        grid on;

        subplot(2, 3, 2);
        plot(frames, root_rot(:,1), 'r-', 'LineWidth', 1.5); hold on;
        plot(frames, root_rot(:,2), 'g-', 'LineWidth', 1.5);
        plot(frames, root_rot(:,3), 'b-', 'LineWidth', 1.5);
        xlabel('Time');
        ylabel('Rotation (degrees)');
        title('Root Rotation Over Time');
        legend('RX', 'RY', 'RZ', 'Location', 'best');
        grid on;

        subplot(2, 3, 3);
        plot(root_pos(:,1), root_pos(:,3), 'k-', 'LineWidth', 1.5); hold on;
        plot(left_ankle_pos(:,1), left_ankle_pos(:,3), 'r-', 'LineWidth', 1.5);
        plot(right_ankle_pos(:,1), right_ankle_pos(:,3), 'b-', 'LineWidth', 1.5);
        plot(left_knee_pos(:,1), left_knee_pos(:,3), 'r:', 'LineWidth', 1.5);
        plot(right_knee_pos(:,1), right_knee_pos(:,3), 'b:', 'LineWidth', 1.5);
        plot(root_pos(1,1), root_pos(1,3), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); 
        plot(root_pos(end,1), root_pos(end,3), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        
        xlabel('Mediolateral (X) Position (m)');
        ylabel('Progression (Z) Position (m)');
        title('Trajectory (Top View / Transverse)');
        legend('Root', 'Left Ankle', 'Right Ankle', 'Left Knee', 'Right Knee', 'Location', 'best');
        grid on;
        axis equal;

        subplot(2, 3, 4);
        plot(root_pos(:,3), root_pos(:,2), 'k-', 'LineWidth', 1.5); hold on;
        plot(left_ankle_pos(:,3), left_ankle_pos(:,2), 'r-', 'LineWidth', 1.5);
        plot(right_ankle_pos(:,3), right_ankle_pos(:,2), 'b-', 'LineWidth', 1.5);
        plot(left_knee_pos(:,3), left_knee_pos(:,2), 'r:', 'LineWidth', 1.5);
        plot(right_knee_pos(:,3), right_knee_pos(:,2), 'b:', 'LineWidth', 1.5);
        for i = 1:inter_space:size(D, 1)
            rx = deg2rad(root_rot(i, 1));
            quiver(root_pos(i,3), root_pos(i,2), cos(rx), sin(rx), 0.1, 'k', 'HandleVisibility', 'off');
        end
        xlabel('Progression (Z) Position (m)');
        ylabel('Height (Y) Position (m)');
        title('Trajectory (Side View / Sagittal)');
        legend('Root', 'Left Ankle', 'Right Ankle', 'Left Knee', 'Right Knee', 'Location', 'best');
        grid on;
        axis equal;

        subplot(2, 3, 5);
        plot(root_pos(:,1), root_pos(:,2), 'k-', 'LineWidth', 1.5); hold on;
        plot(left_ankle_pos(:,1), left_ankle_pos(:,2), 'r-', 'LineWidth', 1.5);
        plot(right_ankle_pos(:,1), right_ankle_pos(:,2), 'b-', 'LineWidth', 1.5);
        plot(left_knee_pos(:,1), left_knee_pos(:,2), 'r:', 'LineWidth', 1.5);
        plot(right_knee_pos(:,1), right_knee_pos(:,2), 'b:', 'LineWidth', 1.5);
        
        xlabel('Mediolateral (X) Position (m)');
        ylabel('Height (Y) Position (m)');
        title('Trajectory (Front View / Frontal)');
        legend('Root', 'Left Ankle', 'Right Ankle', 'Left Knee', 'Right Knee', 'Location', 'best');
        grid on;
        axis equal;

        subplot(2, 3, 6);
        plot3(root_pos(:,1), root_pos(:,3), root_pos(:,2), 'k-', 'LineWidth', 2); hold on;
        plot3(left_ankle_pos(:,1), left_ankle_pos(:,3), left_ankle_pos(:,2), 'r-', 'LineWidth', 2);
        plot3(right_ankle_pos(:,1), right_ankle_pos(:,3), right_ankle_pos(:,2), 'b-', 'LineWidth', 2);
        plot3(left_knee_pos(:,1), left_knee_pos(:,3), left_knee_pos(:,2), 'r:', 'LineWidth', 2);
        plot3(right_knee_pos(:,1), right_knee_pos(:,3), right_knee_pos(:,2), 'b:', 'LineWidth', 2);
        for i = 1:inter_space:size(D, 1)
            u = 0;
            v = cosd(root_rot(i, 1));
            w = sind(root_rot(i, 1));
            quiver3(root_pos(i,1), root_pos(i,3), root_pos(i,2), u, v, w, 0.2, 'k', 'HandleVisibility', 'off');
        end
        xlabel('Mediolateral (X) Position (m)');
        ylabel('Progression (Z) Position (m)');
        zlabel('Height (Y) Position (m)');
        title('3D Trajectory');
        legend('Root', 'Left Ankle', 'Right Ankle', 'Left Knee', 'Right Knee', 'Location', 'best');
        grid on;
        axis equal;
        view(3);

        sgtitle('Root, Knee, and Ankle Motion Analysis', 'FontSize', 16, 'FontWeight', 'bold');

        figure();%2
        plot3(root_pos(:,1), root_pos(:,3), root_pos(:,2), 'k-', 'LineWidth', 2); hold on;
        plot3(left_ankle_pos(:,1), left_ankle_pos(:,3), left_ankle_pos(:,2), 'r-', 'LineWidth', 2);
        plot3(right_ankle_pos(:,1), right_ankle_pos(:,3), right_ankle_pos(:,2), 'b-', 'LineWidth', 2);
        plot3(left_knee_pos(:,1), left_knee_pos(:,3), left_knee_pos(:,2), 'r:', 'LineWidth', 2);
        plot3(right_knee_pos(:,1), right_knee_pos(:,3), right_knee_pos(:,2), 'b:', 'LineWidth', 2);
        plot3(root_pos(1,1), root_pos(1,3), root_pos(1,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
        plot3(root_pos(end,1), root_pos(end,3), root_pos(end,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        
        for i = 1:inter_space:size(D, 1)
            line([left_hip_pos_all(i,1), left_knee_pos(i,1)], [left_hip_pos_all(i,3), left_knee_pos(i,3)], [left_hip_pos_all(i,2), left_knee_pos(i,2)], 'Color', [0.8 0.2 0.2 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            line([left_knee_pos(i,1), left_ankle_pos(i,1)], [left_knee_pos(i,3), left_ankle_pos(i,3)], [left_knee_pos(i,2), left_ankle_pos(i,2)], 'Color', [0.8 0.2 0.2 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            line([right_hip_pos_all(i,1), right_knee_pos(i,1)], [right_hip_pos_all(i,3), right_knee_pos(i,3)], [right_hip_pos_all(i,2), right_knee_pos(i,2)], 'Color', [0.2 0.2 0.8 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            line([right_knee_pos(i,1), right_ankle_pos(i,1)], [right_knee_pos(i,3), right_ankle_pos(i,3)], [right_knee_pos(i,2), right_ankle_pos(i,2)], 'Color', [0.2 0.2 0.8 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            
            u = 0;
            v = cos(left_ankle_orientation(i));
            w = sin(left_ankle_orientation(i));
            quiver3(left_ankle_pos(i,1), left_ankle_pos(i,3), left_ankle_pos(i,2), u, v, w, 0.1, 'r', 'HandleVisibility', 'off');

            u = 0;
            v = cos(right_ankle_orientation(i) );
            w = sin(right_ankle_orientation(i) );
            quiver3(right_ankle_pos(i,1), right_ankle_pos(i,3), right_ankle_pos(i,2), u, v, w, 0.1, 'b', 'HandleVisibility', 'off');

            u = 0;
            v = cosd(root_rot(i, 1));
            w = sind(root_rot(i, 1));
            quiver3(root_pos(i,1), root_pos(i,3), root_pos(i,2), u, v, w, 0.2, 'k', 'HandleVisibility', 'off');
        end

        xlabel('Mediolateral (X) Position (m)');
        ylabel('Progression (Z) Position (m)');
        zlabel('Height (Y) Position (m)');
        title('Detailed 3D Trajectory');
        legend('Root', 'Left Ankle', 'Right Ankle', 'Left Knee', 'Right Knee', 'Start', 'End', 'Location', 'best');
        grid on;
        axis equal;
        view(3);

        figure();%3
        plot(root_pos(:,3), root_pos(:,2), 'k-', 'LineWidth', 1.5); hold on;
        plot(left_ankle_pos(:,3), left_ankle_pos(:,2), 'r-', 'LineWidth', 1.5);
        plot(right_ankle_pos(:,3), right_ankle_pos(:,2), 'b-', 'LineWidth', 1.5);
        plot(left_knee_pos(:,3), left_knee_pos(:,2), 'r:', 'LineWidth', 1.5);
        plot(right_knee_pos(:,3), right_knee_pos(:,2), 'b:', 'LineWidth', 1.5);

        for i = 1:inter_space:size(D, 1)
            line([left_hip_pos_all(i,3), left_knee_pos(i,3)], [left_hip_pos_all(i,2), left_knee_pos(i,2)], 'Color', [0.8 0.2 0.2 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            line([left_knee_pos(i,3), left_ankle_pos(i,3)], [left_knee_pos(i,2), left_ankle_pos(i,2)], 'Color', [0.8 0.2 0.2 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            line([right_hip_pos_all(i,3), right_knee_pos(i,3)], [right_hip_pos_all(i,2), right_knee_pos(i,2)], 'Color', [0.2 0.2 0.8 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            line([right_knee_pos(i,3), right_ankle_pos(i,3)], [right_knee_pos(i,2), right_ankle_pos(i,2)], 'Color', [0.2 0.2 0.8 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');

            u = cos(left_ankle_orientation(i)  );
            v = sin(left_ankle_orientation(i)  );
            quiver(left_ankle_pos(i,3), left_ankle_pos(i,2), u, v, 0.1, 'm', 'HandleVisibility', 'off');

            u = cos(right_ankle_orientation(i)  );
            v = sin(right_ankle_orientation(i)  );
            quiver(right_ankle_pos(i,3), right_ankle_pos(i,2), u, v, 0.1, 'c', 'HandleVisibility', 'off');

            v = cosd(root_rot(i, 1));
            w = sind(root_rot(i, 1));
            quiver(root_pos(i,3), root_pos(i,2), v, w, 0.1, 'm', 'HandleVisibility', 'off');
        end

        xlabel('Progression (Z) Position (m)');
        ylabel('Height (Y) Position (m)');
        title('Ankle and Knee Trajectories (Sagittal View)');
        legend('Root', 'Left Ankle', 'Right Ankle', 'Left Knee', 'Right Knee', 'Location', 'best');
        grid on;
        axis equal;
        
        figure();%4
        hold on;
        plot(left_ankle_pos_FR2(:,3), left_ankle_pos_FR2(:,2), 'r-', 'LineWidth', 1.5);
        plot(right_ankle_pos_FR2(:,3), right_ankle_pos_FR2(:,2), 'b-', 'LineWidth', 1.5);

        for i = 1:inter_space:size(D, 1)
            line([left_hip_pos_all_FR2(i,3), left_knee_pos_FR2(i,3)], [left_hip_pos_all_FR2(i,2), left_knee_pos_FR2(i,2)], 'Color', [0.8 0.2 0.2 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            line([left_knee_pos_FR2(i,3), left_ankle_pos_FR2(i,3)], [left_knee_pos_FR2(i,2), left_ankle_pos_FR2(i,2)], 'Color', [0.8 0.2 0.2 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            line([right_hip_pos_all_FR2(i,3), right_knee_pos_FR2(i,3)], [right_hip_pos_all_FR2(i,2), right_knee_pos_FR2(i,2)], 'Color', [0.2 0.2 0.8 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            line([right_knee_pos_FR2(i,3), right_ankle_pos_FR2(i,3)], [right_knee_pos_FR2(i,2), right_ankle_pos_FR2(i,2)], 'Color', [0.2 0.2 0.8 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');

            u = cos(left_ankle_orientation_FR2(i)  );
            v = sin(left_ankle_orientation_FR2(i)  );
            quiver(left_ankle_pos_FR2(i,3), left_ankle_pos_FR2(i,2), u, v, 0.1, 'm', 'HandleVisibility', 'off');

            u = cos(right_ankle_orientation_FR2(i)  );
            v = sin(right_ankle_orientation_FR2(i)  );
            quiver(right_ankle_pos_FR2(i,3), right_ankle_pos_FR2(i,2), u, v, 0.1, 'c', 'HandleVisibility', 'off');
            
            v = cosd(root_rot(i, 1));
            w = sind(root_rot(i, 1));
            quiver(0,0, v, w, 0.1, 'm', 'HandleVisibility', 'off');
        end

        xlabel('Progression (Z) Position (m)');
        ylabel('Height (Y) Position (m)');
        title('Ankle and Knee Trajectories (Sagittal View)');
        legend('Root', 'Left Ankle', 'Right Ankle', 'Left Knee', 'Right Knee', 'Location', 'best');
        grid on;
        axis equal;
    end

    trajectories.time = frames;
    trajectories.pelvis_orientation = root_rot(:, 1);
    trajectories.left_ankle_pos = left_ankle_pos(:,[3,2]);
    trajectories.right_ankle_pos = right_ankle_pos(:,[3,2]);
    trajectories.left_ankle_pos_FR2 = left_ankle_pos_FR2(:,[3,2]);
    trajectories.right_ankle_pos_FR2 = right_ankle_pos_FR2(:,[3,2]);
    trajectories.left_ankle_orientation = left_ankle_orientation;
    trajectories.right_ankle_orientation = right_ankle_orientation;
    trajectories.left_ankle_orientation_FR2 = left_ankle_orientation_FR2;
    trajectories.right_ankle_orientation_FR2 = right_ankle_orientation_FR2;
end
