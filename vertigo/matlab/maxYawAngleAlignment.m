% Find optimal yaw angle for alignment of the two signals by minimising the
% sum of sqaured differences between points in the two trajectories.
% Note: This assumes the trajectories start at approximately the same time
% (i.e. recording started at approximately the same time).
% INPUT:
% gt : a timeseries object containing time and x,y,z.
% alg: a timeseries object containing the estimate to align to.
% n_yaw_angles : Number of discrete yaw angles to use.
% n_samples    : Number of samples from the algorithm to use.

function yaw = maxYawAngleAlignment(gt, alg, n_yaw_angles, n_samples, plot_figures)

fprintf('Starting yaw angle alignment calculation..\n');

if (nargin < 4)
    plot_figures = false;
end

% Align the start of the data to be at (0,0,0) and t=0.
alg.Data(:,1:3) = alg.Data(:,1:3) - alg.Data(1,1:3);
gt.Data(:,1:3)  = gt.Data(:,1:3)  - gt.Data(1,1:3);
alg.Time = alg.Time - alg.Time(1);
gt.Time  = gt.Time  - gt.Time(1);

% Keep only the first n_samples of the alg to reduce computation time.
alg = delsample(alg,'Index',(n_samples+1):length(alg.Time));
fprintf('Finding optimal yaw angle with %d data points with %d ground truth points\n', length(alg.Time), length(gt.Time));

[alg, gt] = synchronize(alg, gt,'Union','KeepOriginalTimes',true);

fprintf('After synchronization: optimal yaw angle with %d data points with %d ground truth points\n', length(alg.Time), length(gt.Time));

% Find yaw angle that optimally aligns the signals.
initial_guess = pi;
alignment_metric = linspace(0,0,n_yaw_angles);

% For plotting
theta_plot  = [];
metric_plot = [];

n_iterations = 3;

for i=0:(n_iterations-1)
    theta1 = initial_guess - pi/10^i;
    theta2 = initial_guess + pi/10^i;
    theta = linspace(theta1, theta2, n_yaw_angles);
    fprintf('Searching %d angles in the range %f to %f\n', n_yaw_angles, theta1*180/pi, theta2*180/pi);
    
    parfor t = 1:n_yaw_angles
        gt_rot = rotateTs(gt, theta(t));

        % Calculate the metric - for each point, find the distance to the closest point in the other dataset.
        metric = 0;
        for i=1:length(gt_rot.Time)
            diff = alg.Data(:,1:3) - gt_rot.Data(i,1:3);
            dist = diff(:,1).^2 + diff(:,2).^2;
            metric = metric + min(dist);
        end
        alignment_metric(t) = metric;
    end
    
    [~, ind] = min(alignment_metric);
    initial_guess = theta(ind);
    
    theta_plot  = [theta_plot, theta];
    metric_plot = [metric_plot, alignment_metric];
end

[~, ind] = min(alignment_metric);
yaw = theta(ind);

fprintf('Optimal yaw angle: %f degrees\n', yaw*180/pi);

if (plot_figures)
    figure(); subplot(1,2,1); hold on;
    title('Unrotated Data');
    plot(alg.Data(:,1), alg.Data(:,2),'r');
    plot(gt.Data(:,1) , gt.Data(:,2),'k*');
    axis equal
    legend('Algorithm', 'Ground Truth');
    
    subplot(1,2,2); hold on;
    title('Correlation vs theta angle');
    plot(theta_plot*180/pi, metric_plot, '.');
    plot(yaw*180/pi, alignment_metric(ind), '*r');
    xlabel('Rotation around the z axis (deg)');
    ylabel('Sum of distance to closest point (ICP metric)');
end
end