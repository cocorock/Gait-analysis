%% get_amc_files: Finds all AMC files in the 'AMC' directory.
%
% Credits:
%   Victor Ferman, Adrolab FEEC/UNICAMP
%
% Description:
%   This function searches the 'AMC' subdirectory for all files with the .amc
%   extension and returns a structure array of the found files.
%
% Input:
%   None.
%
% Output:
%   amc_files - struct array: A structure array containing information about the AMC files.

function amc_files = get_amc_files()
    fprintf('Searching for AMC files in AMC folder...\n');
    amc_files = dir('AMC/*.amc');
    
    if isempty(amc_files)
        error('No AMC files found in AMC directory!');
    end
    
    fprintf('Found %d AMC files:\n', length(amc_files));
    for i = 1:length(amc_files)
        fprintf('  %d. %s\n', i, amc_files(i).name);
    end
end