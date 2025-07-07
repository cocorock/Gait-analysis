%% Process_Gait_Data_Kinematics_w_Orientation
% This script loads the processed angular data, calculates the linear
% kinematics with orientation, and saves the results.

clear all; close all; clc;

%% Add function paths
addpath('Functions_rev/kinematics_w_orientation');

%% Load Data
fprintf('Loading processed angular data...\n');
load('Gait Data/processed_gait_data_angular_10_samples.mat', 'data');

%% Parameters
phi = 0; % Rotation angle in degrees

%% Calculate Linear Kinematics with Orientation
fprintf('Calculating linear kinematics with orientation...\n');
linear_kinematics_pose = linear_kinematics_w_pose(data, phi);

%% Save Data
fprintf('Saving linear kinematics with orientation data...\n');
save('Gait Data/linear_kinematics_w_orientation.mat', 'linear_kinematics_pose');

fprintf('Processing complete.\n');