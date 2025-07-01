function rotated_points = rotate_data_x90(points)
% Rotates a set of 3D points 90 degrees around the X-axis.
%
% Input:
%   points - An Nx3 matrix where each row is a 3D point (x, y, z).
%
% Output:
%   rotated_points - An Nx3 matrix of the rotated points.

    % Rotation matrix for a +90 degree rotation around the X-axis
    theta = 90;
    R_x = [1, 0, 0;
           0, cosd(theta), -sind(theta);
           0, sind(theta), cosd(theta)];

    % Apply the rotation to each point
    % The points matrix is Nx3, so we transpose it to 3xN for matrix multiplication
    % and then transpose the result back to Nx3.
    rotated_points = (R_x * points')';
end
