function segmented_gait_cycles = segment_gait_cycles(trajectory_struct, left_hs_indices, right_hs_indices)
    % SEGMENT_GAIT_CYCLES Segments a trajectory structure into individual gait cycles.
    %   segmented_gait_cycles = SEGMENT_GAIT_CYCLES(trajectory_struct, left_hs_indices, right_hs_indices)
    %   segments all fields of the input trajectory_struct based on the
    %   detected heel strike indices for left and right legs. A gait cycle
    %   is defined from one heel strike of a limb to the next heel strike
    %   of the same limb.
    %
    %   Inputs:
    %     trajectory_struct: A struct containing various trajectory fields
    %                        (e.g., time, pelvis_orientation, left_ankle_pos, etc.).
    %     left_hs_indices: Indices of left heel strikes.
    %     right_hs_indices: Indices of right heel strikes.
    %
    %   Output:
    %     segmented_gait_cycles: A cell array, where each cell contains a
    %                            struct representing a single gait cycle.

    segmented_gait_cycles = {};
    field_names = fieldnames(trajectory_struct);
    disp(field_names)
    
    % Process left leg gait cycles
    for i = 1:(length(left_hs_indices) - 1)
        start_idx = left_hs_indices(i);
        end_idx = left_hs_indices(i+1);

        current_cycle_struct = struct();
        for j = 1:length(field_names)
            field = field_names{j};
            cycle_data = trajectory_struct.(field)(start_idx:end_idx, :);
            
            if contains(field, 'left', 'IgnoreCase', true) 
                % If the field seems to contain positional data (and not FR2), normalize it
                % by subtracting the last (X,Y) position of the cycle.
                if contains(field, 'pos', 'IgnoreCase', true) && ~contains(field, 'FR2', 'IgnoreCase', true) && size(cycle_data, 2) >= 2
                    last_pos_xy = cycle_data(end, 1:2);
                    cycle_data(:, 1:2) = cycle_data(:, 1:2) - last_pos_xy;
                end
                current_cycle_struct.(field(6:end)) = cycle_data;
                disp('L---')
            end
        end
        
        segmented_gait_cycles{end+1} = current_cycle_struct;
        
    end

    % Process right leg gait cycles
    for i = 1:(length(right_hs_indices) - 1)
        start_idx = right_hs_indices(i);
        end_idx = right_hs_indices(i+1);

        current_cycle_struct = struct();
        for j = 1:length(field_names)
            field = field_names{j};
            cycle_data = trajectory_struct.(field)(start_idx:end_idx, :);
            
           if contains(field, 'right', 'IgnoreCase', true) 
                % If the field seems to contain positional data (and not FR2), normalize it
                % by subtracting the last (X,Y) position of the cycle.
                if contains(field, 'pos', 'IgnoreCase', true) && ~contains(field, 'FR2', 'IgnoreCase', true) && size(cycle_data, 2) >= 2
                    last_pos_xy = cycle_data(end, 1:2);
                    cycle_data(:, 1:2) = cycle_data(:, 1:2) - last_pos_xy;
                end
                current_cycle_struct.(field(7:end)) = cycle_data;
                disp('R---')
            end
        end
        segmented_gait_cycles{end+1} = current_cycle_struct;
        disp(size(segmented_gait_cycles))
    end
end