function rotated_points = rotate_data_xyz(points, angle_x, angle_y, angle_z)
%%
% Rotates a set of 3D points around the X, Y, and Z axes.
%
% Inputs:
%   points    - An Nx3 matrix where each row is a 3D point (x, y, z).
%   angle_x   - Rotation angle around the X-axis in degrees.
%   angle_y   - Rotation angle around the Y-axis in degrees.
%   angle_z   - Rotation angle around the Z-axis in degrees.
%
% Output:
%   rotated_points - An Nx3 matrix of the rotated points.

    % Convert angles to radians
    rad_x = deg2rad(angle_x);
    rad_y = deg2rad(angle_y);
    rad_z = deg2rad(angle_z);

    % Rotation matrices (Z-Y-X Euler angles convention)
    Rz = [cos(rad_z), -sin(rad_z), 0;
          sin(rad_z),  cos(rad_z), 0;
          0,           0,          1];

    Ry = [cos(rad_y),  0, sin(rad_y);
          0,           1, 0;
          -sin(rad_y), 0, cos(rad_y)];

    Rx = [1, 0,          0;
          0, cos(rad_x), -sin(rad_x);
          0, sin(rad_x),  cos(rad_x)];

    % Combined rotation matrix (R = Rx * Ry * Rz)
    R = Rx * Ry * Rz;

    % Apply the rotation to each point
    % The points matrix is Nx3, so we transpose it to 3xN for matrix multiplication
    % and then transpose the result back to Nx3.
    rotated_points = (R * points')';
end