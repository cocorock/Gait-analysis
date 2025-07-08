% This script runs the plot_root_and_ankles_trajectory function
% with the specified ASF and AMC files.
close all;
% Define file paths
asf_file = 'AMC/02.asf';
amc_file = 'AMC/02_01.amc';

% Call the plotting function
plot_root_and_ankles_trajectory(asf_file, amc_file);
