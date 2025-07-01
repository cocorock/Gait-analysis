clear all; clc; close all;
% Load your AMC data first
D = amc_to_matrix('AMC/39_01.amc');

% Plot the root motion
plot_root_motion(D);