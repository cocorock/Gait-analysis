import numpy as np
import scipy.io
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt

def load_gait_data(mat_file_path):
    """
    Loads and preprocesses the complete gait data from the specified .mat file for TP-GMM.
    Includes position, velocity, and orientation for both frames.

    Args:
        mat_file_path (str): The path to the .mat file.

    Returns:
        tuple: (trajectories_fr1, trajectories_fr2, time_data, frame_origins)
               - trajectories_fr1: list of trajectories in frame 1 [time, x, y, vx, vy, orientation]
               - trajectories_fr2: list of trajectories in frame 2 [time, x, y, vx, vy, orientation]  
               - time_data: list of time vectors for each trajectory
               - frame_origins: list of frame origins for each trajectory [(x1,y1), (x2,y2)]
    """
    try:
        mat_data = scipy.io.loadmat(mat_file_path)
        processed_gait_data = mat_data['processed_gait_data']
        
        trajectories_fr1 = []
        trajectories_fr2 = []
        time_data = []
        frame_origins = []
        
        # The data is in a cell array, so we iterate through it
        for i in range(processed_gait_data.shape[1]):
            # The struct is often nested inside a 1x1 array within the cell
            trial_data = processed_gait_data[0, i][0, 0]
            
            # Extract all data fields
            time = trial_data['time'].flatten().astype(np.float64)
            
            # Frame 1 data
            ankle_pos = trial_data['ankle_pos'].astype(np.float64)  # [200x2]
            ankle_vel = trial_data['ankle_pos_velocity'].astype(np.float64)  # [200x2]
            ankle_orient = trial_data['ankle_orientation'].flatten().astype(np.float64)  # [200x1]
            
            # Frame 2 data  
            ankle_pos_fr2 = trial_data['ankle_pos_FR2'].astype(np.float64)  # [200x2]
            ankle_vel_fr2 = trial_data['ankle_pos_FR2_velocity'].astype(np.float64)  # [200x2]
            ankle_orient_fr2 = trial_data['ankle_orientation_FR2'].flatten().astype(np.float64)  # [200x1]
            
            # Estimate frame origins (could be first point, mean, or predefined)
            # For now, we'll use the last point as the origin for each frame
            origin_fr1 = ankle_pos[-1, :]  # Last point of frame 1
            origin_fr2 = ankle_pos_fr2[-1, :]  # Last point of frame 2
            
            # Store frame origins for this trajectory
            frame_origins.append([origin_fr1, origin_fr2])
            
            # Create complete trajectories for both frames
            # Frame 1: [time, x, y, vx, vy, orientation] = 6 dimensions
            trajectory_fr1 = np.column_stack([
                time,
                ankle_pos[:, 0],     # x position
                ankle_pos[:, 1],     # y position  
                ankle_vel[:, 0],     # x velocity
                ankle_vel[:, 1],     # y velocity
                ankle_orient         # orientation
            ])
            
            # Frame 2: [time, x, y, vx, vy, orientation] = 6 dimensions
            trajectory_fr2 = np.column_stack([
                time,
                ankle_pos_fr2[:, 0],     # x position
                ankle_pos_fr2[:, 1],     # y position
                ankle_vel_fr2[:, 0],     # x velocity
                ankle_vel_fr2[:, 1],     # y velocity
                ankle_orient_fr2         # orientation
            ])
            
            trajectories_fr1.append(trajectory_fr1)
            trajectories_fr2.append(trajectory_fr2)
            time_data.append(time)
            
        print(f"Loaded {len(trajectories_fr1)} trajectories")
        print(f"Each trajectory has {trajectory_fr1.shape[1]} dimensions: [time, x, y, vx, vy, orientation]")
        print(f"Trajectory length: {trajectory_fr1.shape[0]} time steps")
        print(f"Frame origins computed for each trajectory")
            
        return trajectories_fr1, trajectories_fr2, time_data, frame_origins
    except FileNotFoundError:
        print(f"Error: The file {mat_file_path} was not found.")
        return None, None, None, None
    except KeyError as e:
        print(f"Error: Could not find required field in the .mat file: {e}")
        return None, None, None, None

class TPGMM:
    """
    Task-Parameterized Gaussian Mixture Model implementation for multi-dimensional gait data.
    Properly handles multiple reference frames and their transformations.
    """
    def __init__(self, n_components=8, n_frames=2):
        self.n_components = n_components
        self.n_frames = n_frames
        self.gmms = []  # One GMM per frame
        self.priors = None
        self.means = None
        self.covariances = None
        self.data_dim = None
        self.frame_origins = None
        
    def fit(self, trajectories_list, frame_origins=None):
        """
        Fit TP-GMM to multiple trajectories from different frames.
        
        Args:
            trajectories_list: List of trajectory lists, one per frame
            frame_origins: Origins of each frame for each trajectory (optional)
        """
        # Store data dimensionality and frame origins
        self.data_dim = trajectories_list[0][0].shape[1]
        self.frame_origins = frame_origins
        print(f"Training TP-GMM with {self.data_dim}D data...")
        
        # Train individual GMMs for each frame
        print("Training individual GMMs for each frame...")
        for frame_idx, trajectories in enumerate(trajectories_list):
            # Concatenate all trajectories for this frame
            frame_data = np.vstack(trajectories)
            print(f"Frame {frame_idx + 1}: {frame_data.shape[0]} data points")
            
            # Train GMM for this frame
            gmm = GaussianMixture(n_components=self.n_components, 
                                covariance_type='full', 
                                random_state=0,
                                max_iter=300,
                                tol=1e-6)
            gmm.fit(frame_data)
            self.gmms.append(gmm)
            print(f"Frame {frame_idx + 1} GMM trained (Log-likelihood: {gmm.score(frame_data):.2f})")
        
        # Store parameters for TP-GMM
        self.priors = self.gmms[0].weights_  # Use first frame's weights as reference
        self.means = [gmm.means_ for gmm in self.gmms]
        self.covariances = [gmm.covariances_ for gmm in self.gmms]
        
        print("TP-GMM training complete.")
        print(f"Model learned {self.n_components} components across {self.n_frames} frames")
    
    def transform_to_frame(self, data, from_frame_origin, to_frame_origin):
        """
        Transform data from one frame to another.
        
        Args:
            data: Data to transform [n_points, n_dims]
            from_frame_origin: Origin of source frame [x, y]
            to_frame_origin: Origin of target frame [x, y]
            
        Returns:
            Transformed data
        """
        # Simple translation transformation (position components only)
        # For more complex transformations, you might need rotation matrices
        transformed_data = data.copy()
        
        # Transform position components (indices 1 and 2 are x,y positions)
        if len(from_frame_origin) >= 2 and len(to_frame_origin) >= 2:
            transformed_data[:, 1] = data[:, 1] - from_frame_origin[0] + to_frame_origin[0]
            transformed_data[:, 2] = data[:, 2] - from_frame_origin[1] + to_frame_origin[1]
        
        return transformed_data
    
    def compute_gaussian_pdf(self, x, mean, cov):
        """
        Compute multivariate Gaussian PDF properly.
        
        Args:
            x: Data point [n_dims]
            mean: Mean vector [n_dims]
            cov: Covariance matrix [n_dims, n_dims]
            
        Returns:
            Gaussian probability density
        """
        n_dims = len(mean)
        
        # Add regularization to avoid singular matrices
        regularized_cov = cov + np.eye(n_dims) * 1e-8
        
        try:
            # Compute inverse and determinant
            cov_inv = np.linalg.inv(regularized_cov)
            cov_det = np.linalg.det(regularized_cov)
            
            # Handle negative determinant
            if cov_det <= 0:
                return 1e-10
            
            # Compute difference
            diff = x - mean
            
            # Compute PDF
            exp_term = np.exp(-0.5 * diff.T @ cov_inv @ diff)
            norm_term = 1.0 / np.sqrt((2 * np.pi) ** n_dims * cov_det)
            
            return norm_term * exp_term
            
        except np.linalg.LinAlgError:
            return 1e-10
    
    def predict_proba_tp(self, input_data, frame_weights, frame_origins=None):
        """
        Predict component probabilities for TP-GMM using proper product of Gaussians.
        
        Args:
            input_data: Input data points [n_points, input_dim]
            frame_weights: Weights for each frame [n_frames]
            frame_origins: Origins for each frame (optional)
            
        Returns:
            Component probabilities [n_points, n_components]
        """
        n_points = input_data.shape[0]
        input_dim = input_data.shape[1]
        
        # Initialize probabilities
        component_probs = np.zeros((n_points, self.n_components))
        
        print(f"Computing TP-GMM probabilities for {n_points} points...")
        
        for i in range(n_points):
            point = input_data[i]
            
            # Compute responsibility for each component
            for k in range(self.n_components):
                
                # PRODUCT OF GAUSSIANS across frames
                # This is the key TP-GMM operation
                gaussian_product = 1.0
                
                for frame_idx in range(self.n_frames):
                    if frame_weights[frame_idx] > 0:
                        # Get parameters for this frame and component
                        mean_k = self.means[frame_idx][k]
                        cov_k = self.covariances[frame_idx][k]
                        
                        # Extract input dimensions
                        mu_I = mean_k[:input_dim]
                        Sigma_II = cov_k[:input_dim, :input_dim]
                        
                        # Compute Gaussian PDF for this frame
                        gaussian_val = self.compute_gaussian_pdf(point, mu_I, Sigma_II)
                        
                        # Raise to power of frame weight (for product of Gaussians)
                        # This is the proper TP-GMM formulation
                        gaussian_product *= gaussian_val ** frame_weights[frame_idx]
                
                # Multiply by prior probability
                component_probs[i, k] = self.priors[k] * gaussian_product
        
        # Normalize probabilities for each point
        for i in range(n_points):
            total = np.sum(component_probs[i, :])
            if total > 0:
                component_probs[i, :] /= total
            else:
                component_probs[i, :] = 1.0 / self.n_components
        
        print("TP-GMM probabilities computed.")
        return component_probs
    
    def perform_gmr(self, input_data, frame_weights, target_frame_origin=None):
        """
        Perform Gaussian Mixture Regression with proper task parameterization.
        
        Args:
            input_data: Input data [n_points, input_dim] (typically time)
            frame_weights: Weights for each frame [n_frames]
            target_frame_origin: Origin for the target frame (optional)
            
        Returns:
            Predicted output data [n_points, output_dim]
        """
        input_dim = input_data.shape[1]
        output_dim = self.data_dim - input_dim
        
        print(f"Performing TP-GMR with frame weights: {frame_weights}")
        
        # Get component probabilities using proper product of Gaussians
        component_probs = self.predict_proba_tp(input_data, frame_weights)
        
        # Perform GMR
        expected_output = np.zeros((input_data.shape[0], output_dim))
        
        for i in range(input_data.shape[0]):
            input_point = input_data[i]
            gmr_point = np.zeros(output_dim)
            
            for k in range(self.n_components):
                # PRODUCT OF CONDITIONAL DISTRIBUTIONS
                # This is the proper TP-GMM conditional expectation
                
                # Compute weighted parameters across frames
                weighted_mean = np.zeros(output_dim)
                weighted_precision = np.zeros((output_dim, output_dim))
                total_weight = 0
                
                for frame_idx in range(self.n_frames):
                    if frame_weights[frame_idx] > 0:
                        # Get parameters for this frame and component
                        mean_k = self.means[frame_idx][k]
                        cov_k = self.covariances[frame_idx][k]
                        
                        # Partition parameters
                        mu_I = mean_k[:input_dim]
                        mu_O = mean_k[input_dim:]
                        Sigma_II = cov_k[:input_dim, :input_dim]
                        Sigma_OO = cov_k[input_dim:, input_dim:]
                        Sigma_OI = cov_k[input_dim:, :input_dim]
                        
                        # Add regularization
                        Sigma_II += np.eye(input_dim) * 1e-8
                        Sigma_OO += np.eye(output_dim) * 1e-8
                        
                        try:
                            # Conditional parameters
                            Sigma_II_inv = np.linalg.inv(Sigma_II)
                            conditional_mean = mu_O + Sigma_OI @ Sigma_II_inv @ (input_point - mu_I)
                            conditional_cov = Sigma_OO - Sigma_OI @ Sigma_II_inv @ Sigma_OI.T
                            
                            # Product of Gaussians: combine precisions
                            conditional_precision = np.linalg.inv(conditional_cov + np.eye(output_dim) * 1e-8)
                            
                            # Weight by frame importance
                            weight = frame_weights[frame_idx]
                            weighted_mean += weight * conditional_precision @ conditional_mean
                            weighted_precision += weight * conditional_precision
                            total_weight += weight
                            
                        except np.linalg.LinAlgError:
                            # Fallback to simple mean
                            weighted_mean += frame_weights[frame_idx] * mu_O
                            total_weight += frame_weights[frame_idx]
                
                # Compute final conditional expectation
                if total_weight > 0:
                    try:
                        if np.linalg.det(weighted_precision) > 1e-10:
                            final_mean = np.linalg.inv(weighted_precision) @ weighted_mean
                        else:
                            final_mean = weighted_mean / total_weight
                    except np.linalg.LinAlgError:
                        final_mean = weighted_mean / total_weight
                else:
                    final_mean = np.zeros(output_dim)
                
                # Add to GMR result weighted by component probability
                gmr_point += component_probs[i, k] * final_mean
            
            expected_output[i] = gmr_point
        
        print("TP-GMR completed.")
        return expected_output

def plot_comprehensive_results(trajectories_fr1, trajectories_fr2, tpgmm, time_vector):
    """
    Create comprehensive plots showing position, velocity, and orientation results.
    """
    frame_combinations = [
        ([1.0, 0.0], "Frame 1 Only"),
        ([0.0, 1.0], "Frame 2 Only"), 
        ([0.5, 0.5], "Equal Frames"),
        ([0.7, 0.3], "Mostly Frame 1"),
        ([0.3, 0.7], "Mostly Frame 2")
    ]
    
    fig, axes = plt.subplots(3, len(frame_combinations), figsize=(20, 12))
    
    for col, (frame_weights, label) in enumerate(frame_combinations):
        print(f"Generating plots for {label}...")
        reproduced_trajectory = tpgmm.perform_gmr(time_vector, frame_weights)
        
        # Position plot (top row)
        ax = axes[0, col]
        # Plot original trajectories (faded)
        for traj in trajectories_fr1:
            ax.plot(traj[:, 1], traj[:, 2], 'b-', alpha=0.1, linewidth=0.5)
        for traj in trajectories_fr2:
            ax.plot(traj[:, 1], traj[:, 2], 'g-', alpha=0.1, linewidth=0.5)
        
        # Plot reproduced trajectory
        ax.plot(reproduced_trajectory[:, 0], reproduced_trajectory[:, 1], 
                'r-', linewidth=3, label='TP-GMR')
        ax.set_title(f'Position: {label}')
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.grid(True)
        ax.axis('equal')
        
        # Velocity plot (middle row)
        ax = axes[1, col]
        # Plot velocity magnitude over time
        time_flat = time_vector.flatten()
        vel_magnitude = np.sqrt(reproduced_trajectory[:, 2]**2 + reproduced_trajectory[:, 3]**2)
        ax.plot(time_flat, vel_magnitude, 'r-', linewidth=3, label='TP-GMR Velocity')
        
        # Plot original velocity magnitudes for comparison
        for traj in trajectories_fr1:
            vel_mag_orig = np.sqrt(traj[:, 3]**2 + traj[:, 4]**2)
            ax.plot(traj[:, 0], vel_mag_orig, 'b-', alpha=0.1, linewidth=0.5)
        for traj in trajectories_fr2:
            vel_mag_orig = np.sqrt(traj[:, 3]**2 + traj[:, 4]**2)
            ax.plot(traj[:, 0], vel_mag_orig, 'g-', alpha=0.1, linewidth=0.5)
            
        ax.set_title(f'Velocity Magnitude: {label}')
        ax.set_xlabel('Time')
        ax.set_ylabel('Velocity Magnitude')
        ax.grid(True)
        
        # Orientation plot (bottom row)
        ax = axes[2, col]
        ax.plot(time_flat, reproduced_trajectory[:, 4], 'r-', linewidth=3, label='TP-GMR Orientation')
        
        # Plot original orientations for comparison
        for traj in trajectories_fr1:
            ax.plot(traj[:, 0], traj[:, 5], 'b-', alpha=0.1, linewidth=0.5)
        for traj in trajectories_fr2:
            ax.plot(traj[:, 0], traj[:, 5], 'g-', alpha=0.1, linewidth=0.5)
            
        ax.set_title(f'Orientation: {label}')
        ax.set_xlabel('Time')
        ax.set_ylabel('Orientation (rad)')
        ax.grid(True)
    
    plt.tight_layout()
    plt.show()

def main():
    """
    Main function to run the proper TP-GMM/GMR example.
    """
    # --- 1. Load and Prepare Data ---
    mat_file = 'Gait Data/new_processed_gait_data.mat'
    trajectories_fr1, trajectories_fr2, time_data, frame_origins = load_gait_data(mat_file)

    if trajectories_fr1 is None:
        return

    # --- 2. Train the Task-Parameterized Gaussian Mixture Model ---
    tpgmm = TPGMM(n_components=6, n_frames=2)  # Reduced components for stability
    tpgmm.fit([trajectories_fr1, trajectories_fr2], frame_origins)

    # --- 3. Perform Task-Parameterized GMR ---
    n_repro_points = 200
    time_vector = np.linspace(0, 1, n_repro_points)[:, np.newaxis]

    # Create comprehensive plots
    plot_comprehensive_results(trajectories_fr1, trajectories_fr2, tpgmm, time_vector)

    # --- 4. Demonstrate Product of Gaussians ---
    print("\n=== Product of Gaussians Demonstration ===")
    
    # Test with different frame weights to show the effect
    test_weights = [
        [1.0, 0.0],
        [0.0, 1.0], 
        [0.5, 0.5],
        [0.8, 0.2],
        [0.2, 0.8]
    ]
    
    # Sample a few time points
    test_times = np.array([[0.2], [0.5], [0.8]])
    
    for weights in test_weights:
        probs = tpgmm.predict_proba_tp(test_times, weights)
        print(f"Frame weights {weights}: Component probabilities at t=[0.2,0.5,0.8]")
        print(f"  Component distribution: {np.mean(probs, axis=0)}")
        print()

    # --- 5. Frame Origin Analysis ---
    print("=== Frame Origins Analysis ===")
    if frame_origins:
        print(f"Number of trajectories with frame origins: {len(frame_origins)}")
        
        # Analyze the difference between frame origins
        origin_differences = []
        for origins in frame_origins:
            diff = np.linalg.norm(origins[1] - origins[0])
            origin_differences.append(diff)
        
        print(f"Average distance between frame origins: {np.mean(origin_differences):.3f}")
        print(f"Min/Max distance between frame origins: {np.min(origin_differences):.3f}/{np.max(origin_differences):.3f}")
        
        # Plot frame origins
        plt.figure(figsize=(10, 6))
        origins_fr1 = np.array([origins[0] for origins in frame_origins])
        origins_fr2 = np.array([origins[1] for origins in frame_origins])
        
        plt.scatter(origins_fr1[:, 0], origins_fr1[:, 1], c='blue', label='Frame 1 Origins', alpha=0.6)
        plt.scatter(origins_fr2[:, 0], origins_fr2[:, 1], c='green', label='Frame 2 Origins', alpha=0.6)
        
        # Draw lines connecting corresponding origins
        for i in range(len(frame_origins)):
            plt.plot([origins_fr1[i, 0], origins_fr2[i, 0]], 
                    [origins_fr1[i, 1], origins_fr2[i, 1]], 
                    'k-', alpha=0.3, linewidth=0.5)
        
        plt.xlabel('X Position')
        plt.ylabel('Y Position')
        plt.title('Frame Origins Comparison')
        plt.legend()
        plt.grid(True)
        plt.axis('equal')
        plt.show()
    
    print("\nComplete TP-GMM analysis finished!")

if __name__ == '__main__':
    main()