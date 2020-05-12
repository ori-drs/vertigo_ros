% This function aligns one trajectory to another in yaw and time.
% INPUT:
% alg: timeseries containing trajectory to align ground truth to
% gt : timeseries containing ground truth to be aligned.
% OUTPUT:
% gt_aligned: timeseries with aligned ground truth data.

function gt_aligned = alignData(alg, gt, plot_figs)

if (nargin<3)
    plot_figs = false;
end

fprintf('Aligning %d data points with %d ground truth points\n', length(alg.Time), length(gt.Time));

% Variables
n_states_for_yaw_alignment  = 400*20; % Number of alg states used for yaw alignment.
n_states_for_time_alignment = 100;   % Number of alg states used for time alignment.
n_yaw_angles                = 100;    % Number of discrete yaw angles.

% Yaw alignment
yaw = maxYawAngleAlignment(gt, alg, n_yaw_angles, n_states_for_yaw_alignment, plot_figs);

gt_rot = rotateTs(gt, yaw);
gt_rot.Data(:,1:3) = gt_rot.Data(:,1:3) + alg.Data(1,1:3);

if (plot_figs)
    figure(); subplot(1,2,1); hold on;
    title('Unrotated Data');
    plot(alg.Data(:,1) - alg.Data(1,1), alg.Data(:,2) - alg.Data(1,2),'r');
    plot(gt.Data(:,1)  - gt.Data(1,1) , gt.Data(:,2) - gt.Data(1,2)  ,'k');
    axis equal
    legend('Algorithm', 'Ground Truth');
    
    subplot(1,2,2); hold on; title('Rotated trajectories');
    plot(alg.Data(:,1), alg.Data(:,2),'r');
    plot(gt_rot.Data(:,1) , gt_rot.Data(:,2),'k');
    axis equal
    legend('Algorithm', 'Ground Truth');
    
    figure(); hold on;
    title('Rotated trajectories in 3D');
    plot3(alg.Data(:,1), alg.Data(:,2), alg.Data(:,3),'r');
    plot3(gt_rot.Data(:,1), gt_rot.Data(:,2), gt_rot.Data(:,3),'k');
    axis equal
    legend('Algorithm', 'Ground Truth');    
end

% Temporal alignment
%t_offset = maxTimeAlignment(gt_rot, alg, n_states_for_time_alignment, true);
t_offset = 0;

gt_aligned = gt_rot;
gt_aligned.Time = gt_aligned.Time - t_offset;

if (plot_figs)
    figure(); title('Final yaw and time-aligned data');
    label = 'xyz';
    for i=1:3
        ax(i) = subplot(3,1,i); hold on; title(label(i));
        plot(gt_aligned.Time, gt_aligned.Data(:,i),'k.');
        plot(alg.Time, alg.Data(:,i),'r');
    end
    linkaxes([ax(1) ax(2) ax(3)], 'xy');
end

end