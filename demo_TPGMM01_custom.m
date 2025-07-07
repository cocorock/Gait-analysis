function demo_TPGMM01_custom
% Custom Task-parameterized Gaussian mixture model (TP-GMM) encoding.
% Adapted for gait analysis data with position and velocity.

% addpath('./m_fcts/'); % Assuming m_fcts is in the same directory as the original demo
addpath('./Gait Data/'); % To load your TPGMM_data.mat

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.nbStates = 9; %Number of Gaussians in the GMM
model.nbFrames = 2; %Number of candidate frames of reference
model.nbVar = 4; %Dimension of the datapoints in the dataset (here: x,y,vx,vy)


%% Load 3rd order tensor data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load custom 3rd order tensor data...');
load('TPGMM_data.mat'); % Load your 's' structure and 'nbSamples'

% Regenerate data tensor for TP-GMM learning
nbData = s(1).nbData; % Should be 200
Data = zeros(model.nbVar, model.nbFrames, nbSamples*nbData);

for n=1:nbSamples
    % s(n).Data contains [pos_x; pos_y; vel_x; vel_y]
    % We need to transform both position and velocity components.
    % The transformation for position is (A \ (pos - b))
    % The transformation for velocity is (A \ vel)
    
    for m=1:model.nbFrames
        % Extract transformation components for this frame
        A_inv = inv(s(n).p(m).A); % Inverse of rotation matrix
        b = s(n).p(m).b;          % Translation vector
        
        % Original position and velocity data for this demonstration
        pos_orig = s(n).Data(1:2,:);
        vel_orig = s(n).Data(3:4,:);
        
        % Transform position: Rotate and then translate
        pos_transformed = A_inv * (pos_orig - repmat(b, 1, nbData));
        
        % Transform velocity: Only rotate (velocity vectors are not affected by translation)
        vel_transformed = A_inv * vel_orig;
        
        % Combine transformed position and velocity
        Data(:,m,(n-1)*nbData+1:n*nbData) = [pos_transformed; vel_transformed];
    end
end


%% TP-GMM learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameters estimation of TP-GMM with EM:\n');
model = init_tensorGMM_kmeans(Data, model); 
model = EM_tensorGMM(Data, model);

% Reconstruct GMM for each demonstration
for n=1:nbSamples
	[s(n).Mu, s(n).Sigma] = productTPGMM0(model, s(n).p);
end


%% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position',[10,10,2300,900]);
xx = round(linspace(1,64,nbSamples));
clrmap = colormap('jet');
clrmap = min(clrmap(xx,:),.95);
limAxes = [-1.2 0.8 -1.1 0.9]; % Adjust as needed for your data range


% DEMOS (Original data in world frame)
subplot(1,model.nbFrames+1,1); hold on; box on; title('Demonstrations (Original)');
for n=1:nbSamples
	%Plot frames

	%Plot trajectories (only position components)
	plot(s(n).Data(1,1), s(n).Data(2,1),'.','markersize',15,'color',clrmap(n,:));
	plot(s(n).Data(1,:), s(n).Data(2,:),'-','linewidth',1.5,'color',clrmap(n,:));
	%Plot Gaussians (these are from the product of transformed Gaussians, so they are in the world frame)
	plotGMM(s(n).Mu(1:2,:), s(n).Sigma(1:2,1:2,:), [.5 .5 .5],.8); % Only plot position part of GMM
end
axis(limAxes); axis square; set(gca,'xtick',[],'ytick',[]);

% MODEL (Transformed data in each frame)
p0.A = eye(2);
p0.b = zeros(2,1);
for m=1:model.nbFrames
	subplot(1,model.nbFrames+1,1+m); hold on; grid on; box on; title(['Frame ' num2str(m) ' (Transformed)']);
	for n=1:nbSamples
		% Plot transformed trajectories (only position components)
		plot(squeeze(Data(1,m,(n-1)*nbData+1)), ...
			squeeze(Data(2,m,(n-1)*nbData+1)), '.','markersize',15,'color',clrmap(n,:));
		plot(squeeze(Data(1,m,(n-1)*nbData+1:n*nbData)), ...
			squeeze(Data(2,m,(n-1)*nbData+1:n*nbData)), '-','linewidth',1.5,'color',clrmap(n,:));
	end
	% Plot GMMs for each frame (only position components)
	plotGMM(squeeze(model.Mu(1:2,m,:)), squeeze(model.Sigma(1:2,1:2,m,:)), [.5 .5 .5],.8);

	axis equal; axis([-4.5 4.5 -1 8]); set(gca,'xtick',[0],'ytick',[0]); % Adjust as needed
end

pause;
close all;
end

% Function to plot pegs (copied from original demo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = plotPegs(p, colPegs, fa)
	if ~exist('colPegs')
		colPegs = [0.2863 0.0392 0.2392; 0.9137 0.4980 0.0078];
	end
	if ~exist('fa')
		fa = .6;
	end
	pegMesh = [-4 -3.5; -4 10; -1.5 10; -1.5 -1; 1.5 -1; 1.5 10; 4 10; 4 -3.5; -4 -3.5]' *1E-1;
	for m=1:length(p)
		dispMesh = p(m).A * pegMesh + repmat(p(m).b,1,size(pegMesh,2));
		h(m) = patch(dispMesh(1,:),dispMesh(2,:),colPegs(m,:),'linewidth',1,'edgecolor','none','facealpha',fa);
	end
end