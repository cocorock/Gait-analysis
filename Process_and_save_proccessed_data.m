%% Process and Save All AMC Files
% This script performs an enhanced gait analysis on AMC files. It processes
% multiple AMC files, extracts gait cycles, filters the data

close all; % Close all existing figures
clc;       % Clear the command window
clear; % Clear all variables from the workspace

% Add necessary paths to access functions and data
addpath('./AMC/');
addpath('./Functions_rev/'); % Use the revised functions
addpath('./Gait Data/');

%% Create Output Directories
% This function creates the necessary output directories for saving plots and data.
create_output_directorie();

%% Main Execution
fprintf('=== Gait Analysis v5 Started ===\n');

% Get all AMC files from the specified directory.
amc_files = get_amc_files();

% Process all AMC files to extract and collect gait cycle data.
[all_cycles_data, file_info] = process_all_amc_files(amc_files, true, true);

% Apply filtering to the collected gait data and calculate derivatives.
processed_data = apply_filtering_and_derivatives(all_cycles_data);

% Save the processed data in multiple formats for further analysis.
save_processed_data(processed_data, file_info, 10);

% Reshape the processed data for further analysis and plotting.
reshaped_data = reshapeProcessedData(processed_data);

% Plot the joint cycles from the reshaped data.
plot_joint_cycles(reshaped_data);

fprintf('\n===  Process and Save ===\n');
