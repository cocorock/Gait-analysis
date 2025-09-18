% Setup Environment
% Clear workspace, close figures, and add necessary paths.

close all;
clc;
clear;

% Add necessary paths to access functions and data
addpath('./AMC/');
addpath('./Functions_rev/'); % Use the revised functions
addpath('./Gait Data/');
addpath('./Functions_rev/Version2_4D/')
%% Create Output Directories
% This function creates the necessary output directories for saving plots and data

create_output_directorie();

%% Main Execution - Data Processing (V3)
% Start the gait analysis process.

fprintf('\n=== Process Plot and Save All AMC Files Completed V3 ===\n');

% Get all AMC files from the specified directory.
amc_files = get_amc_files();
% Process all AMC files to extract and collect gait cycle data.
[all_cycles_data, file_info] = process_all_amc_files_v2(amc_files, false, true);

% --- MODIFICATION: Use only cycles from the right leg criterion ---
% Clear the left leg data to ensure all subsequent analysis uses only the
% cycles segmented based on the right leg's heel strike.
all_cycles_data.left_leg_cycles = [];
all_cycles_data.file_indices.left = [];
fprintf('Modified: Using only %d cycles from the right leg criterion.\n', length(all_cycles_data.right_leg_cycles));
% ----------------------------------------------------------------

% Apply filtering to the collected gait data (derivatives are NOT calculated here).
processed_data = apply_filtering_V3(all_cycles_data);

% Ploting right leg angular kinematics (positions only).
plot_gait_kinematics_v3(processed_data, 'right');

% % Ploting left leg angular kinematics
%  plot_gait_kinematics_v3(processed_data, 'left');

% Save the processed data in multiple formats for further analysis.
% NOTE: This is commented out as processed_data in V3 only contains filtered
% positions, which may not be what the saving function expects.
% save_processed_data_4D(processed_data);

%% Kinematics Calculation and Plotting (V3)
% Calculate and visualize linear kinematics.

% Calculate linear kinematics using the new function.
% This function now calculates angular derivatives internally.
linear_kinematics = calculate_linear_kinematics_v3(processed_data, -90);

save_linear_kinematics_structuredV3(linear_kinematics, processed_data.time_standard, file_info);

%  %Plotting Kinematics
plot_linear_kinematics_positionsV2(linear_kinematics, processed_data.time_standard);
% 
% plot_linear_kinematics_velocities(linear_kinematics, processed_data.time_standard);
% 
% plot_linear_kinematics_accelerations(linear_kinematics, processed_data.time_standard);

plot_linear_kinematics_positions_xy(linear_kinematics);

plot_linear_kinematics_velocities_xy(linear_kinematics);

plot_linear_kinematics_accelerations_xy(linear_kinematics);
