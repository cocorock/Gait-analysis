clear all; clc; close all;
% Load your AMC data first
D = amc_to_matrix('39_14.amc');

% Plot the root motion
plot_root_motion(D);