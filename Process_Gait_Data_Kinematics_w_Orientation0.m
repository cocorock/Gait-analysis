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
linear_kinematics_pose = linear_kinematics_pose(1:6);
%% Save Data
% fprintf('Saving linear kinematics with orientation data...\n');
% save('Gait Data/linear_kinematics_w_orientation.mat', 'linear_kinematics_pose');
% fprintf('Processing complete.\n');

%% Transforming
frame_ref_index = 200;
[demos200, frames200] = create_TPGMM_data_FR1(linear_kinematics_pose, frame_ref_index);
frame_ref_index = 120;
[demos120, frames120] = create_TPGMM_data_FR1(linear_kinematics_pose, frame_ref_index);

[demosST, framesST] = create_TPGMM_data_FR_start(linear_kinematics_pose);

[demos_relative_start, frames_relative_start] = create_TPGMM_data_FR_relative_start(linear_kinematics_pose);

%% Plotting section
plot_TPGMM_data(demos200, frames200, linear_kinematics_pose);

plot_TPGMM_data(demos120, frames120, linear_kinematics_pose);

plot_TPGMM_data(demos_relative_start, frames_relative_start, linear_kinematics_pose);

plot_kinematics_with_orientation(linear_kinematics_pose);


%% Save Data for TP-GMM
% fprintf('\nSaving TP-GMM compatible data s to Gait Data/TPGMM_data.mat...\n');
% save('Gait Data/TPGMM_data.mat', 's', 'nbSamples');
% fprintf('TP-GMM data processing complete.\n');




