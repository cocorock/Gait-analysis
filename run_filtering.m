try
    addpath('./filtering');
    create_filtering_plot;
catch ME
    fprintf('Error running script: %s\n', ME.message);
    exit(1);
end
exit(0);