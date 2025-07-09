% This script runs the plot_root_and_ankles_trajectory function
% with the specified ASF and AMC files.
close all;

% Add necessary paths to access functions and data
addpath('./AMC/');
addpath('./Functions_rev/'); % Use the revised functions

% Define file paths
asf_file = 'AMC/02.asf';
amc_file = 'AMC/02_01.amc';

% Call the plotting function
tra = plot_root_and_ankles_trajectory(asf_file, amc_file, true);


% R_lfemur = eul2rotm(deg2rad([90,0,0]), 'XYZ')                      

% extract_hip_knee_flexion(amc_file);
