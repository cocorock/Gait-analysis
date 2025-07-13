
% Functions_rev/calculate_velocity.m
function velocity = calculate_velocity(position_data, dt)
    % CALCULATE_VELOCITY Calculates velocity from position data.
    %   velocity = CALCULATE_VELOCITY(position_data, dt) calculates the
    %   velocity using the formula [pos(i+1)-pos(i)]/dt. Circular
    %   referencing is used for the first and last elements.
    %
    %   Inputs:
    %     position_data: NxM matrix of position data, where N is the number
    %                    of samples and M is the number of dimensions.
    %     dt: Time step between samples (1/sampling_frequency).
    %
    %   Output:
    %     velocity: NxM matrix of calculated velocities.

    if nargin < 2
        error('Not enough input arguments. Usage: calculate_velocity(position_data, dt)');
    end

    num_samples = size(position_data, 1);
    num_dimensions = size(position_data, 2);
    velocity = zeros(num_samples, num_dimensions);
%     fprintf('velocity %d,%d  dt: %2.4f\n', num_samples, num_dimensions, dt);
    for i = 1:num_samples
        if i == 1
            % Forward difference for first element
            velocity(i, :) = (position_data(i+1, :) - position_data(i, :)) / dt;
        elseif i == num_samples
            % Backward difference for last element
            velocity(i, :) = (position_data(i, :) - position_data(i-1, :)) / dt;
        else
            % Central difference for interior points
            velocity(i, :) = (position_data(i+1, :) - position_data(i-1, :)) / (2*dt);
        end
    end
end
