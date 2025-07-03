import scipy.io
import numpy as np
import pandas as pd
import os

# Input .mat file
mat_file_path = 'Gait Data/linear_kinematics_gait_cycles_10_samples.mat'

# Load the .mat file
mat_data = scipy.io.loadmat(mat_file_path)

# The data is stored in a variable named 'data', which is a cell array.
# In Python, this becomes a numpy array of objects, where each object is a demo.
data_cell = mat_data['data']

# Get the number of demonstrations (gait cycles)
n_demos = data_cell.shape[1]

# Assuming the number of points is the number of columns in the first demo
# and all demos have the same number of points.
n_points = data_cell[0, 0].shape[1]

# Create a time vector, similar to the example
t = np.linspace(0, 2 * np.pi, n_points)

trajectories = []
for i in range(n_demos):
    # Get the data for the i-th demo
    demo_data = data_cell[0, i]

    # Extract position and velocity
    pos_x = demo_data[0, :]
    pos_y = demo_data[1, :]
    vel_x = demo_data[2, :]
    vel_y = demo_data[3, :]

    # Combine into a trajectory array: [time, pos_x, pos_y, vel_x, vel_y]
    # The shape should be (n_points, 5)
    traj = np.column_stack([t, pos_x, pos_y, vel_x, vel_y])
    trajectories.append(traj)

# Combine all trajectories into a single numpy array, similar to 'all_data'
all_data = np.vstack(trajectories)

# You can now use 'all_data' which has a shape of (n_demos * n_points, 5)
print("Shape of the final data:", all_data.shape)

# --- Save to CSV with a similar name in the same directory ---
# Get the directory and the filename from the input path
input_dir, input_filename = os.path.split(mat_file_path)
# Get the filename without the extension
base_filename = os.path.splitext(input_filename)[0]
# Create the output filename
output_filename = f"{base_filename}.csv"
# Create the full output path
output_csv_path = os.path.join(input_dir, output_filename)

# Save to the new CSV file path
np.savetxt(output_csv_path, all_data, delimiter=',')
print(f"Data saved to {output_csv_path}")

# To use with pandas:
df = pd.DataFrame(all_data, columns=['time', 'pos_x', 'pos_y', 'vel_x', 'vel_y'])
print("Pandas DataFrame head:")
print(df.head())
