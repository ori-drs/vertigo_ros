% Find optimal time shift for alignment to move the GT signal to the
% estimated signal time.
% Note: This function ASSUMES the signals are spatially aligned.
% INPUT:
% gt: ground truth time series
% alg: algorithm estimate time series
% N: discretization of time shift.

function t_offset = maxTimeAlignment(gt, alg, N, plot_figs)

if (nargin < 4)
    plot_figs = false;
end

% Align first 33% of data before the values start to diverge too much.
percent = 0.33;
gt_range  = floor(length(gt.Time)*percent) :length(gt.Time);
alg_range = floor(length(alg.Time)*percent):length(alg.Time);
gt  = delsample(gt ,'Index', gt_range);
alg = delsample(alg,'Index', alg_range);

fprintf('Starting time alignment calculation with %f oercent of the data..\n', percent*100);

% Find the minimum and maximum time offsets where the signals still have
% non-zero overlap.
min_dt = gt.Time(1) - alg.Time(end);
max_dt = gt.Time(end) - alg.Time(1);

MIN_OVERLAP = 0.8;

dt_plot = [];
corr_plot = [];

n_iterations = 3;

for i=1:n_iterations
    fprintf('Searching %d samples between %f and %f seconds\n', N, min_dt, max_dt);

    dt   = linspace(min_dt, max_dt, N);
    corr = linspace(0,0,N);

    parfor i=1:N 
        % Find the indices where the GT signal overlaps with the estimated
        % signal, taking into account the proposed time delay.
        gt_time_shifted = gt.Time - dt(i);
        overlap_indices = find( (gt_time_shifted > alg.Time(1)) .* ...
                                (gt_time_shifted < alg.Time(end)));
        if ( length(overlap_indices) / length(gt_time_shifted) < MIN_OVERLAP )
            continue
        end

        gt_shifted = timeseries(gt.Data(overlap_indices,:), gt_time_shifted(overlap_indices));

        % Interpolate to match the lower frequency ground truth signal.
        [t, ind] = unique(alg.Time);
        x = zeros(length(gt_shifted.Time), 3);
        x(:,1) = interp1(t, alg.Data(ind,1), gt_shifted.Time );
        x(:,2) = interp1(t, alg.Data(ind,2), gt_shifted.Time );
        x(:,3) = interp1(t, alg.Data(ind,3), gt_shifted.Time );

        % Find the correlation of the two signals.
        % Only correlate x and y as they drift less than z.
        corr(i) = norm( (gt_shifted.Data(:,1) - x(:,1)).^2 + ...
                        (gt_shifted.Data(:,2) - x(:,2)).^2);  
        % corr(i) = norm( (gt_shifted.Data(:,1) - x(:,1)).^2 + ...
        %                 (gt_shifted.Data(:,2) - x(:,2)).^2 + ...
        %                 (gt_shifted.Data(:,3) - x(:,3)).^2 );
    end
    
    non_zero_indices = find(corr > 0);
    dt = dt(non_zero_indices);
    corr = corr(non_zero_indices);
    [~, ind] = min(corr);
    t_offset = dt(ind);
    min_dt = t_offset - (max_dt - min_dt)/10^i;
    max_dt = t_offset + (max_dt - min_dt)/10^i;
    fprintf('Best time offset so far: %f seconds\n', t_offset);
    
    dt_plot   = [dt_plot, dt];
    corr_plot = [corr_plot, corr];
end

fprintf('Optimal time offset: %f seconds\n', t_offset);

if (plot_figs)
    figure(); hold on;
    title('Correlation vs timeshift');
    plot(dt_plot, corr_plot, '.');
    plot(t_offset, corr(ind), 'r*');
    xlabel('delta t (sec)');
    ylabel('correlation');
end