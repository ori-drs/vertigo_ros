% Script to align the trajectories from the GT and data.
% INPUT:
% slam result.csvin tum format
% ground truth.csv in tum format
% output filename: name of output csv file.
%
% Note: All computation is done in the Leica Prism frame (due to the
% ambiguity in transforming from the prism to the base frame).

clc; close all; 

data_filename = 'slam_result_tum.csv';
gt_filename  = 'gt_tum.csv';
output_filename = 'slam_result_tum_aligned.csv';
plot_figs       = true;

% The guard enable you to run the script multiple times without reloading
% data.
fprintf('Loading gt data..\n');
if (~exist('gt_ts'))
    gt_ts = loadData(gt_filename,'gt')
end
fprintf('Loaded %d GT data points\n', length(gt_ts.Time));

if (~exist('data_ts'))
    data_ts = loadData(data_filename,'data')
end
fprintf('Loaded %d data points\n', length(data_ts.Time));

data_aligned_ts = alignData(gt_ts, data_ts, plot_figs);

% Save to file.
t = data_aligned_ts.Time;
d = data_aligned_ts.Data;
csv_output = zeros(8, length(d));
csv_output(1,:) = t;
csv_output(2,:) = d(:,1);
csv_output(3,:) = d(:,2);
csv_output(4,:) = d(:,3);
% Add identity quaternion as this is unknown.
csv_output(5,:) = d(:,5);
csv_output(6,:) = d(:,6);
csv_output(7,:) = d(:,7);
csv_output(8,:) = d(:,4);

csv_output = array2table(csv_output');

writetable(csv_output, output_filename,'Delimiter',' ');