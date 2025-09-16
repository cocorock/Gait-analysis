% close all; % Close all existing figures
clc;       % Clear the command window
clear; % Clear all variables from the workspace
close all;

% Add necessary paths to access functions and data
addpath('./Functions_rev/'); % Use the revised functions
addpath('./Gait Data/');
addpath('./AMC2/');

% Get all AMC files from the specified directory.
amc_files = get_amc_files('AMC2/');

% Initialize a cell array to store trajectories for each file
all_trajectories = cell(length(amc_files), 1);

% Loop through each amc file
for i = 1:length(amc_files)
    % Get the current amc file name
    amc_file = amc_files(i).name;
    
    % Construct the asf file name from the amc file name
    asf_file = [amc_file(1:2) '.asf'];
    
    fprintf('Processing file %d of %d: %s with %s\n', i, length(amc_files), amc_file, asf_file);
    
    % Get the trajectories without plotting the figures
    trajectories = plot_root_and_ankles_trajectoryFIXED(asf_file, amc_file, true, true, false);
    
    % Store the trajectories in the cell array
    all_trajectories{i} = trajectories;
end

% flag10 =  false;
% for i=1:size(amc_files,1)
%     if contains(amc_files(i).name, '10.amc', 'IgnoreCase', true) 
%         flag10 = true;
%     end
% end
flag10 =  true;

  if ~flag10  
    output_file = sprintf('./Gait Data/all_trajectories_ALL-10#%s.mat', amc_file(1:2));
    save(output_file, 'all_trajectories');
  else
    output_file = sprintf('./Gait Data/all_trajectories_ALL#%s.mat', amc_file(1:2));
    save(output_file, 'all_trajectories');
  end

fprintf('All trajectories saved to %s\n', output_file);

