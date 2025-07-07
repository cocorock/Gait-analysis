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

%% Create Data for Task-Parameterized GMM (TP-GMM)
fprintf('Creating data structure for TP-GMM...\n');

num_cycles = length(linear_kinematics_pose);
demos = cell(1, num_cycles);
% Frames will store the transformation for each demonstration's coordinate system
% We will store the origin position of each new frame.
% A common way to store frames is using homogeneous transformation matrices.
% For a simple translation, the matrix is [I, P; 0, 1], where P is the translation vector.
frames = zeros(3, 3, num_cycles); 

for i = 1:num_cycles
    % Get the data for the current cycle
    cycle_data = linear_kinematics_pose{i};
    
    % The new frame of reference is at the end oftrans the trajectory
    origin_pos = cycle_data.pos(:, end);
    
    % Store the transformation matrix for this frame
    % This represents the position of the new frame's origin w.r.t the world frame
    frames(:, :, i) = [eye(2), origin_pos; 0 0 1];
    
    % Transform the position data to be relative to the new origin
    transformed_pos = cycle_data.pos - origin_pos;
    
    % Combine all data for the demonstration.
    % Velocity, acceleration, and orientation are properties of the motion
    % itself and do not change with a simple translation of the coordinate system.
    demos{i} = [transformed_pos; 
                cycle_data.vel; 
                cycle_data.acc;
                cycle_data.orientation;
                cycle_data.orientation_vel;
                cycle_data.orientation_acc];
end

%% Save Data for TP-GMM
fprintf('Saving data for TP-GMM...\n');
save('Gait Data/TPGMM_data.mat', 'demos', 'frames');

fprintf('TP-GMM data processing complete.\n');

%% Plot the TP-GMM data
plot_TPGMM_data(demos, frames, linear_kinematics_pose);

%% Plot the results for the first gait cycle
plot_kinematics_with_orientation(linear_kinematics_pose);