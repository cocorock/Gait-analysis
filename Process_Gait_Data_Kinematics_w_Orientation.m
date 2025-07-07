%% Process_Gait_Data_Kinematics_w_Orientation
% This script loads the processed angular data, calculates the linear
% kinematics with orientation, and saves the results.

clear all; close all; clc;

%% Add function paths
addpath('Functions_rev/kinematics_w_orientation');
addpath('./Functions_rev/'); 
addpath('./Gait Data/');

%% Load Data
fprintf('Loading processed angular data...\n');
load('Gait Data/processed_gait_data_fields_angular.mat');

Data_structured = create_demos_structure_per_cycleV2(mat_data.processed_data);

%% Parameters
phi = -90; % Rotation angle in degrees
%% Calculate Linear Kinematics with Orientation
fprintf('Calculating linear kinematics with orientation...\n');
linear_kinematics_pose = linear_kinematics_w_pose(Data_structured, phi);

%% Save Data
fprintf('Saving linear kinematics with orientation data...\n');
save('Gait Data/linear_kinematics_w_orientation.mat', 'linear_kinematics_pose');

fprintf('Processing complete.\n');

%% Create Data for Task-Parameterized GMM (TP-GMM) compatible with demo_TPGMM01
fprintf('\nCreating data structure s for TP-GMM demo compatibility...\n');

nbSamples = 5; % Use only the first 5 gait cycles
nbData = 200; % Number of data points per cycle
s = struct(); % Initialize empty struct array

for n = 1:nbSamples
    cycle_data = linear_kinematics_pose{n};
    
    % s(n).Data: [pos_x; pos_y; vel_x; vel_y]
    s(n).Data = [cycle_data.pos; cycle_data.vel];
    s(n).nbData = nbData;
    
    % --- Define Frame 1 (index 200) ---
    frame_index_1 = 200;
    origin_pos_1 = cycle_data.pos(:, frame_index_1);
    final_orientation_1 = cycle_data.orientation(frame_index_1);
    R1 = [cos(final_orientation_1), -sin(final_orientation_1);
          sin(final_orientation_1),  cos(final_orientation_1)];
    s(n).p(1).b = origin_pos_1; % Translation vector
    s(n).p(1).A = R1;           % Rotation matrix
    
    % --- Define Frame 2 (index 160) ---
    frame_index_2 = 160;
    origin_pos_2 = cycle_data.pos(:, frame_index_2);
    final_orientation_2 = cycle_data.orientation(frame_index_2);
    R2 = [cos(final_orientation_2), -sin(final_orientation_2);
          sin(final_orientation_2),  cos(final_orientation_2)];
    s(n).p(2).b = origin_pos_2;
    s(n).p(2).A = R2;
end

%% Save Data for TP-GMM
fprintf('\nSaving TP-GMM compatible data s to Gait Data/TPGMM_data.mat...\n');
save('Gait Data/TPGMM_data.mat', 's', 'nbSamples');

fprintf('TP-GMM data processing complete.\n');

%% Plotting section
plot_TPGMM_s_data(s, nbSamples);

plot_kinematics_with_orientation(linear_kinematics_pose);

% plot_TPGMM_data(demos, frames, linear_kinematics_pose);