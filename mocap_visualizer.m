function mocap_visualizer(asf_file, amc_file)
    % MOCAP_VISUALIZER - Reads ASF/AMC files and visualizes walking motion
    % Usage: mocap_visualizer('skeleton.asf', 'walking.amc')
    
    if nargin < 2
        error('Please provide both ASF and AMC file paths');
    end
    
    % Read ASF file (skeleton definition)
    fprintf('Reading ASF file: %s\n', asf_file);
    skeleton = read_asf(asf_file);
    
    % Read AMC file (motion data)
    fprintf('Reading AMC file: %s\n', amc_file);
    motion = read_amc(amc_file, skeleton);
    
    % Extract lower limb sizes
    lower_limb_sizes = extract_lower_limb_sizes(skeleton);
    fprintf('Lower limb sizes extracted:\n');
    fprintf('Hip width: %.2f\n', lower_limb_sizes.hip_width);
    fprintf('Left femur length: %.2f\n', lower_limb_sizes.left_femur);
    fprintf('Right femur length: %.2f\n', lower_limb_sizes.right_femur);
    fprintf('Left tibia length: %.2f\n', lower_limb_sizes.left_tibia);
    fprintf('Right tibia length: %.2f\n', lower_limb_sizes.right_tibia);
    
    % Calculate joint positions for all frames
    fprintf('Calculating joint positions...\n');
    joint_positions = calculate_joint_positions(skeleton, motion);
    
    % Create 3D visualization
    fprintf('Creating 3D visualization...\n');
    visualize_walking_motion(joint_positions, skeleton);
    
    % Create 2D plane projections
    fprintf('Creating 2D plane projections...\n');
    plot_orthogonal_planes(joint_positions);
end

function skeleton = read_asf(filename)
    % Read ASF file and extract skeleton structure
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open ASF file: %s', filename);
    end
    
    skeleton = struct();
    skeleton.joints = containers.Map();
    skeleton.joint_order = {};
    
    current_section = '';
    
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)
            line = strtrim(line);
            
            % Skip empty lines and comments
            if isempty(line) || line(1) == '#'
                continue;
            end
            
            % Check for section headers
            if line(1) == ':'
                current_section = line(2:end);
                continue;
            end
            
            % Parse based on current section
            switch current_section
                case 'units'
                    if contains(line, 'mass')
                        skeleton.mass_unit = extract_unit_value(line);
                    elseif contains(line, 'length')
                        skeleton.length_unit = extract_unit_value(line);
                    elseif contains(line, 'angle')
                        skeleton.angle_unit = extract_unit_value(line);
                    end
                    
                case 'root'
                    tokens = strsplit(line);
                    if strcmp(tokens{1}, 'position')
                        skeleton.root_position = [str2double(tokens{2}), ...
                                                str2double(tokens{3}), ...
                                                str2double(tokens{4})];
                    elseif strcmp(tokens{1}, 'orientation')
                        skeleton.root_orientation = [str2double(tokens{2}), ...
                                                   str2double(tokens{3}), ...
                                                   str2double(tokens{4})];
                    end
                    
                case 'bonedata'
                    if strcmp(line, 'begin')
                        joint = struct();
                        joint_name = '';
                    elseif strcmp(line, 'end')
                        if ~isempty(joint_name)
                            skeleton.joints(joint_name) = joint;
                            skeleton.joint_order{end+1} = joint_name;
                        end
                    else
                        tokens = strsplit(line);
                        if length(tokens) >= 2
                            switch tokens{1}
                                case 'id'
                                    joint.id = str2double(tokens{2});
                                case 'name'
                                    joint_name = tokens{2};
                                    joint.name = joint_name;
                                case 'direction'
                                    joint.direction = [str2double(tokens{2}), ...
                                                     str2double(tokens{3}), ...
                                                     str2double(tokens{4})];
                                case 'length'
                                    joint.length = str2double(tokens{2});
                                case 'axis'
                                    joint.axis = [str2double(tokens{2}), ...
                                                str2double(tokens{3}), ...
                                                str2double(tokens{4})];
                                case 'dof'
                                    joint.dof = tokens(2:end);
                                case 'limits'
                                    % Parse limits for each DOF
                                    limits = [];
                                    for i = 2:2:length(tokens)
                                        if i+1 <= length(tokens)
                                            limits = [limits; str2double(tokens{i}(2:end)), ...
                                                    str2double(tokens{i+1}(1:end-1))];
                                        end
                                    end
                                    joint.limits = limits;
                            end
                        end
                    end
                    
                case 'hierarchy'
                    if strcmp(line, 'begin')
                        skeleton.hierarchy = containers.Map();
                    elseif ~strcmp(line, 'end')
                        tokens = strsplit(line);
                        parent = tokens{1};
                        children = tokens(2:end);
                        skeleton.hierarchy(parent) = children;
                    end
            end
        end
    end
    
    fclose(fid);
end

function unit_value = extract_unit_value(line)
    % Extract unit value from a line like "mass 1.0"
    tokens = strsplit(line);
    if length(tokens) >= 2
        unit_value = str2double(tokens{2});
    else
        unit_value = 1.0;
    end
end

function motion = read_amc(filename, skeleton)
    % Read AMC file and extract motion data
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open AMC file: %s', filename);
    end
    
    motion = struct();
    motion.frames = {};
    frame_count = 0;
    
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)
            line = strtrim(line);
            
            % Skip empty lines and comments
            if isempty(line) || line(1) == '#' || line(1) == ':'
                continue;
            end
            
            % Check if line is a frame number
            if ~isnan(str2double(line))
                frame_count = frame_count + 1;
                motion.frames{frame_count} = struct();
                motion.frames{frame_count}.number = str2double(line);
                continue;
            end
            
            % Parse joint data
            if frame_count > 0
                tokens = strsplit(line);
                joint_name = tokens{1};
                
                if skeleton.joints.isKey(joint_name) || strcmp(joint_name, 'root')
                    values = [];
                    for i = 2:length(tokens)
                        values = [values, str2double(tokens{i})];
                    end
                    motion.frames{frame_count}.(joint_name) = values;
                end
            end
        end
    end
    
    motion.num_frames = frame_count;
    fclose(fid);
end

function lower_limb_sizes = extract_lower_limb_sizes(skeleton)
    % Extract sizes of lower limb segments based on VICON structure
    
    lower_limb_sizes = struct();
    
    % Initialize with default values
    lower_limb_sizes.hip_width = 0;
    lower_limb_sizes.left_femur = 0;
    lower_limb_sizes.right_femur = 0;
    lower_limb_sizes.left_tibia = 0;
    lower_limb_sizes.right_tibia = 0;
    
    % Length scale factor from ASF units
    length_scale = 1/0.45*25.4/1000; % Based on your ASF units (0.45 scale factor)
    
    % Extract femur lengths using exact joint names from your ASF
    if skeleton.joints.isKey('lfemur')
        joint = skeleton.joints('lfemur');
        if isfield(joint, 'length')
            lower_limb_sizes.left_femur = joint.length * length_scale;
        end
    end
    
    if skeleton.joints.isKey('rfemur')
        joint = skeleton.joints('rfemur');
        if isfield(joint, 'length')
            lower_limb_sizes.right_femur = joint.length * length_scale;
        end
    end
    
    % Extract tibia lengths
    if skeleton.joints.isKey('ltibia')
        joint = skeleton.joints('ltibia');
        if isfield(joint, 'length')
            lower_limb_sizes.left_tibia = joint.length * length_scale;
        end
    end
    
    if skeleton.joints.isKey('rtibia')
        joint = skeleton.joints('rtibia');
        if isfield(joint, 'length')
            lower_limb_sizes.right_tibia = joint.length * length_scale;
        end
    end
    
    % Calculate hip width from hip joint positions
    if skeleton.joints.isKey('lhipjoint') && skeleton.joints.isKey('rhipjoint')
        lhip_joint = skeleton.joints('lhipjoint');
        rhip_joint = skeleton.joints('rhipjoint');
        
        if isfield(lhip_joint, 'length') && isfield(rhip_joint, 'length')
            % Hip width is approximately the sum of both hip joint lengths
            lower_limb_sizes.hip_width = (lhip_joint.length + rhip_joint.length) * length_scale;
        end
    end
    
    % If hip width is still zero, estimate from directions
    if lower_limb_sizes.hip_width == 0
        lower_limb_sizes.hip_width = 200; % Default estimate in mm
    end
end

function joint_positions = calculate_joint_positions(skeleton, motion)
    % Calculate 3D positions of joints for all frames
    
    num_frames = motion.num_frames;
    joint_positions = struct();
    
    % Initialize position arrays for key joints
    joint_positions.root = zeros(num_frames, 3);
    joint_positions.lhip = zeros(num_frames, 3);
    joint_positions.rhip = zeros(num_frames, 3);
    joint_positions.lknee = zeros(num_frames, 3);
    joint_positions.rknee = zeros(num_frames, 3);
    joint_positions.lankle = zeros(num_frames, 3);
    joint_positions.rankle = zeros(num_frames, 3);
    
    for frame = 1:num_frames
        frame_data = motion.frames{frame};
        
        % Start with root position and orientation
        root_pos = [0, 0, 0];
        root_rot = [0, 0, 0];
        
        if isfield(frame_data, 'root')
            root_data = frame_data.root;
            if length(root_data) >= 6
                root_pos = root_data(1:3);
                root_rot = root_data(4:6);
            end
        end
        
        joint_positions.root(frame, :) = root_pos;
        
        % Calculate positions using forward kinematics
        joint_positions = forward_kinematics(skeleton, frame_data, joint_positions, frame, root_pos, root_rot);
    end
end

function joint_positions = forward_kinematics(skeleton, frame_data, joint_positions, frame, root_pos, root_rot)
    % Perform forward kinematics to calculate joint positions using VICON joint names
    
    % Create root transformation matrix
    T_root = create_transform_matrix(root_pos, root_rot);
    
    % Calculate hip positions using hip joint data
    hip_offset = 100; % mm - approximate hip width
    
    % Get hip joint lengths and directions if available
    lhip_offset = [-hip_offset/2, 0, 0];
    rhip_offset = [hip_offset/2, 0, 0];
    
    if skeleton.joints.isKey('lhipjoint')
        lhip_joint = skeleton.joints('lhipjoint');
        if isfield(lhip_joint, 'direction') && isfield(lhip_joint, 'length')
            lhip_offset = lhip_joint.direction * lhip_joint.length * 45; % Apply scale
        end
    end
    
    if skeleton.joints.isKey('rhipjoint')
        rhip_joint = skeleton.joints('rhipjoint');
        if isfield(rhip_joint, 'direction') && isfield(rhip_joint, 'length')
            rhip_offset = rhip_joint.direction * rhip_joint.length * 45; % Apply scale
        end
    end
    
    % Transform hip offsets by root transformation
    lhip_world = T_root * [lhip_offset'; 1];
    rhip_world = T_root * [rhip_offset'; 1];
    
    joint_positions.lhip(frame, :) = lhip_world(1:3)';
    joint_positions.rhip(frame, :) = rhip_world(1:3)';
    
    % Calculate left leg chain
    joint_positions = calculate_leg_chain(skeleton, frame_data, joint_positions, frame, ...
                                        joint_positions.lhip(frame, :), 'left');
    
    % Calculate right leg chain
    joint_positions = calculate_leg_chain(skeleton, frame_data, joint_positions, frame, ...
                                        joint_positions.rhip(frame, :), 'right');
end

function joint_positions = calculate_leg_chain(skeleton, frame_data, joint_positions, frame, hip_pos, side)
    % Calculate knee and ankle positions for one leg using proper joint names
    
    % Default segment lengths (in mm, scaled by length unit)
    length_scale = 1/0.45*25.4/1000; % Based on your ASF units (0.45 scale factor)
    default_femur = 40 * length_scale;
    default_tibia = 35 * length_scale;
    
    % Get actual segment lengths from skeleton
    femur_length = default_femur;
    tibia_length = default_tibia;
    
    % Joint names based on your ASF structure
    if strcmp(side, 'left')
        femur_joint_name = 'lfemur';
        tibia_joint_name = 'ltibia';
        foot_joint_name = 'lfoot';
    else
        femur_joint_name = 'rfemur';
        tibia_joint_name = 'rtibia';
        foot_joint_name = 'rfoot';
    end
    
    % Get segment lengths from skeleton
    if skeleton.joints.isKey(femur_joint_name)
        joint = skeleton.joints(femur_joint_name);
        if isfield(joint, 'length') && joint.length > 0
            femur_length = joint.length * length_scale;
        end
    end
    
    if skeleton.joints.isKey(tibia_joint_name)
        joint = skeleton.joints(tibia_joint_name);
        if isfield(joint, 'length') && joint.length > 0
            tibia_length = joint.length * length_scale;
        end
    end
    
    % Get joint angles from frame data
    femur_angles = [0, 0, 0];
    tibia_angles = [0, 0, 0];
    foot_angles = [0, 0, 0];
    
    if isfield(frame_data, femur_joint_name)
        femur_data = frame_data.(femur_joint_name);
        if length(femur_data) >= 3
            femur_angles = femur_data(1:3);
        elseif length(femur_data) == 1
            femur_angles(1) = femur_data(1);
        end
    end
    
    if isfield(frame_data, tibia_joint_name)
        tibia_data = frame_data.(tibia_joint_name);
        if length(tibia_data) >= 1
            tibia_angles(1) = tibia_data(1); % Tibia typically has only one DOF (knee flexion)
        end
    end
    
    if isfield(frame_data, foot_joint_name)
        foot_data = frame_data.(foot_joint_name);
        if length(foot_data) >= 2
            foot_angles(1:2) = foot_data(1:2);
        elseif length(foot_data) == 1
            foot_angles(1) = foot_data(1);
        end
    end
    
    % Get joint directions and axes from skeleton
    femur_direction = [0, -1, 0]; % Default downward
    tibia_direction = [0, -1, 0]; % Default downward
    
    if skeleton.joints.isKey(femur_joint_name)
        joint = skeleton.joints(femur_joint_name);
        if isfield(joint, 'direction')
            femur_direction = joint.direction;
        end
    end
    
    if skeleton.joints.isKey(tibia_joint_name)
        joint = skeleton.joints(tibia_joint_name);
        if isfield(joint, 'direction')
            tibia_direction = joint.direction;
        end
    end
    
    % Calculate knee position
    femur_transform = create_transform_matrix([0,0,0], femur_angles);
    % Create homogeneous vector: direction * length as 3D vector, then add 1 for homogeneous coordinate
    femur_vec_3d = femur_direction * femur_length;
    femur_vector_homo = [femur_vec_3d(:); 1];  % Ensure column vector
    femur_vector = femur_transform * femur_vector_homo;
    knee_pos = hip_pos + femur_vector(1:3)';
    
    % Calculate ankle position
    tibia_transform = create_transform_matrix([0,0,0], tibia_angles);
    combined_transform = femur_transform * tibia_transform;
    % Create homogeneous vector for tibia
    tibia_vec_3d = tibia_direction * tibia_length;
    tibia_vector_homo = [tibia_vec_3d(:); 1];  % Ensure column vector
    tibia_vector = combined_transform * tibia_vector_homo;
    ankle_pos = hip_pos + tibia_vector(1:3)';
    
    % Store positions
    if strcmp(side, 'left')
        joint_positions.lknee(frame, :) = knee_pos;
        joint_positions.lankle(frame, :) = ankle_pos;
    else
        joint_positions.rknee(frame, :) = knee_pos;
        joint_positions.rankle(frame, :) = ankle_pos;
    end
end

function T = create_transform_matrix(translation, rotation)
    % Create 4x4 transformation matrix from translation and rotation
    
    % Convert angles to radians
    rx = deg2rad(rotation(1));
    ry = deg2rad(rotation(2));
    rz = deg2rad(rotation(3));
    
    % Create rotation matrices
    Rx = [1, 0, 0; 0, cos(rx), -sin(rx); 0, sin(rx), cos(rx)];
    Ry = [cos(ry), 0, sin(ry); 0, 1, 0; -sin(ry), 0, cos(ry)];
    Rz = [cos(rz), -sin(rz), 0; sin(rz), cos(rz), 0; 0, 0, 1];
    
    % Combined rotation matrix
    R = Rz * Ry * Rx;
    
    % Create transformation matrix
    T = eye(4);
    T(1:3, 1:3) = R;
    T(1:3, 4) = translation;
end

function visualize_walking_motion(joint_positions, skeleton)
    % Create 3D visualization of walking motion
    
    figure('Name', 'Walking Motion Visualization', 'Position', [100, 100, 1200, 800]);
    
    % Create subplot for 3D animation
    subplot(2, 2, [1, 3]);
    
    num_frames = size(joint_positions.root, 1);
    
    % Plot trajectories
    hold on;
    plot3(joint_positions.root(:,1), joint_positions.root(:,2), joint_positions.root(:,3), ...
          'b-', 'LineWidth', 2, 'DisplayName', 'Root');
    plot3(joint_positions.lknee(:,1), joint_positions.lknee(:,2), joint_positions.lknee(:,3), ...
          'r-', 'LineWidth', 1.5, 'DisplayName', 'Left Knee');
    plot3(joint_positions.rknee(:,1), joint_positions.rknee(:,2), joint_positions.rknee(:,3), ...
          'g-', 'LineWidth', 1.5, 'DisplayName', 'Right Knee');
    plot3(joint_positions.lankle(:,1), joint_positions.lankle(:,2), joint_positions.lankle(:,3), ...
          'm-', 'LineWidth', 1.5, 'DisplayName', 'Left Ankle');
    plot3(joint_positions.rankle(:,1), joint_positions.rankle(:,2), joint_positions.rankle(:,3), ...
          'c-', 'LineWidth', 1.5, 'DisplayName', 'Right Ankle');
    
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Z (mm)');
    title('3D Walking Motion Trajectories');
    legend('Location', 'best');
    grid on;
    axis equal;
    view(45, 30);
    
    % Create additional plots for analysis
    subplot(2, 2, 2);
    plot(1:num_frames, joint_positions.root(:,3), 'b-', 'LineWidth', 2);
    hold on;
    plot(1:num_frames, joint_positions.lknee(:,3), 'r-', 'LineWidth', 1.5);
    plot(1:num_frames, joint_positions.rknee(:,3), 'g-', 'LineWidth', 1.5);
    xlabel('Frame');
    ylabel('Height (Z) mm');
    title('Vertical Motion');
    legend('Root', 'Left Knee', 'Right Knee', 'Location', 'best');
    grid on;
    
    subplot(2, 2, 4);
    plot(1:num_frames, joint_positions.root(:,1), 'b-', 'LineWidth', 2);
    hold on;
    plot(1:num_frames, joint_positions.lankle(:,1), 'm-', 'LineWidth', 1.5);
    plot(1:num_frames, joint_positions.rankle(:,1), 'c-', 'LineWidth', 1.5);
    xlabel('Frame');
    ylabel('Forward Distance (X) mm');
    title('Forward Progression');
    legend('Root', 'Left Ankle', 'Right Ankle', 'Location', 'best');
    grid on;
    
    hold off;
end

function plot_orthogonal_planes(joint_positions)
    % Plot motion trajectories in X-Y, Y-Z, and X-Z planes
    
    figure('Name', 'Orthogonal Plane Projections', 'Position', [150, 150, 1200, 800]);
    
    % Define colors for each joint
    colors = struct();
    colors.root = [0, 0, 1];        % Blue
    colors.lhip = [1, 0.5, 0];      % Orange
    colors.rhip = [0.5, 0.5, 0.5];  % Gray
    colors.lknee = [1, 0, 0];       % Red
    colors.rknee = [0, 0.5, 0];     % Dark Green
    colors.lankle = [1, 0, 1];      % Magenta
    colors.rankle = [0, 1, 1];      % Cyan
    
    % Joint names for legend
    joint_names = {'root', 'lhip', 'rhip', 'lknee', 'rknee', 'lankle', 'rankle'};
    legend_names = {'Root', 'Left Hip', 'Right Hip', 'Left Knee', 'Right Knee', 'Left Ankle', 'Right Ankle'};
    
    % X-Y Plane (Top view - transverse plane)
    subplot(2, 2, 1);
    hold on;
    for i = 1:length(joint_names)
        joint_name = joint_names{i};
        if isfield(joint_positions, joint_name)
            plot(joint_positions.(joint_name)(:,1), joint_positions.(joint_name)(:,2), ...
                 'Color', colors.(joint_name), 'LineWidth', 2, 'DisplayName', legend_names{i});
            % Add start and end markers
            plot(joint_positions.(joint_name)(1,1), joint_positions.(joint_name)(1,2), ...
                 'o', 'Color', colors.(joint_name), 'MarkerSize', 8, 'MarkerFaceColor', colors.(joint_name));
            plot(joint_positions.(joint_name)(end,1), joint_positions.(joint_name)(end,2), ...
                 's', 'Color', colors.(joint_name), 'MarkerSize', 8, 'MarkerFaceColor', colors.(joint_name));
        end
    end
    xlabel('X (mm) - Forward/Backward');
    ylabel('Y (mm) - Left/Right');
    title('X-Y Plane (Top View - Transverse)');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    axis equal;
    
    % Y-Z Plane (Front view - coronal plane)
    subplot(2, 2, 2);
    hold on;
    for i = 1:length(joint_names)
        joint_name = joint_names{i};
        if isfield(joint_positions, joint_name)
            plot(joint_positions.(joint_name)(:,2), joint_positions.(joint_name)(:,3), ...
                 'Color', colors.(joint_name), 'LineWidth', 2, 'DisplayName', legend_names{i});
            % Add start and end markers
            plot(joint_positions.(joint_name)(1,2), joint_positions.(joint_name)(1,3), ...
                 'o', 'Color', colors.(joint_name), 'MarkerSize', 8, 'MarkerFaceColor', colors.(joint_name));
            plot(joint_positions.(joint_name)(end,2), joint_positions.(joint_name)(end,3), ...
                 's', 'Color', colors.(joint_name), 'MarkerSize', 8, 'MarkerFaceColor', colors.(joint_name));
        end
    end
    xlabel('Y (mm) - Left/Right');
    ylabel('Z (mm) - Up/Down');
    title('Y-Z Plane (Front View - Coronal)');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    axis equal;
    
    % X-Z Plane (Side view - sagittal plane)
    subplot(2, 2, 3);
    hold on;
    for i = 1:length(joint_names)
        joint_name = joint_names{i};
        if isfield(joint_positions, joint_name)
            plot(joint_positions.(joint_name)(:,1), joint_positions.(joint_name)(:,3), ...
                 'Color', colors.(joint_name), 'LineWidth', 2, 'DisplayName', legend_names{i});
            % Add start and end markers
            plot(joint_positions.(joint_name)(1,1), joint_positions.(joint_name)(1,3), ...
                 'o', 'Color', colors.(joint_name), 'MarkerSize', 8, 'MarkerFaceColor', colors.(joint_name));
            plot(joint_positions.(joint_name)(end,1), joint_positions.(joint_name)(end,3), ...
                 's', 'Color', colors.(joint_name), 'MarkerSize', 8, 'MarkerFaceColor', colors.(joint_name));
        end
    end
    xlabel('X (mm) - Forward/Backward');
    ylabel('Z (mm) - Up/Down');
    title('X-Z Plane (Side View - Sagittal)');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    axis equal;
    
    % Combined trajectory plot with different line styles
    subplot(2, 2, 4);
    hold on;
    
    % Plot only key joints with different line styles for clarity
    key_joints = {'root', 'lknee', 'rknee', 'lankle', 'rankle'};
    key_names = {'Root', 'Left Knee', 'Right Knee', 'Left Ankle', 'Right Ankle'};
    line_styles = {'-', '--', '-.', ':', '-'};
    
    for i = 1:length(key_joints)
        joint_name = key_joints{i};
        if isfield(joint_positions, joint_name)
            % Plot X coordinate vs frame number
            plot(1:size(joint_positions.(joint_name), 1), joint_positions.(joint_name)(:,1), ...
                 line_styles{i}, 'Color', colors.(joint_name), 'LineWidth', 2, ...
                 'DisplayName', [key_names{i} ' (X)']);
        end
    end
    
    xlabel('Frame Number');
    ylabel('X Position (mm)');
    title('Forward Progression Over Time');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    
    hold off;
end