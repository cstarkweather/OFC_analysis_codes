% plot PSTH for approach and avoid trials
apr_x = [];apr_y = [];
av_x = [];av_y = [];
for i = 1:11

    for n = 1:n_trials

        if trial_type(i).decisions(n)==1% & trial_type(i).decision_times(n)>1500
            decision_time = round(trial_type(i).decision_times(n));
            apr_x = [apr_x; zeros(1,6000-decision_time) trajectory(i).trajectories_x(n,1:decision_time)];
            apr_y = [apr_y; zeros(1,6000-decision_time) trajectory(i).trajectories_y(n,1:decision_time)];
        else
            if trial_type(i).decisions(n)==0 %& trial_type(i).decision_times(n)>1500
            decision_time = round(trial_type(i).decision_times(n));
            av_x = [av_x; zeros(1,6000-decision_time) trajectory(i).trajectories_x(n,1:decision_time)];
            av_y = [av_y; zeros(1,6000-decision_time) trajectory(i).trajectories_y(n,1:decision_time)];
            end
            end

    end


end
%%
summed_apr_x = sum(apr_x);summed_apr_y = sum(apr_y);
for i=1:6000
    denom_apr(i) = size(apr_x,1) - sum(apr_x(:,i)==0);
end



summed_av_x = sum(av_x);summed_av_y = sum(av_y);
for i=1:6000
    denom_av(i) = size(av_x,1) - sum(av_x(:,i)==0);
end

smoothing = 1;
subplot(1,2,1)
plot(smooth(summed_apr_x./denom_apr,smoothing),'b')
hold on
plot(smooth(summed_av_x./denom_av,smoothing),'r')
ylim([0 5])
xlim([2000 6000])

subplot(1,2,2)
plot(smooth(summed_apr_y./denom_apr,smoothing),'b--')
hold on
plot(smooth(summed_av_y./denom_av,smoothing),'r--')
ylim([0 5])
xlim([2000 6000])
%%
% Example script to plot PSTHs with SEM shading.

% ------------------ Parameters --------------------------------
smoothing = 100;                  % Your smoothing window
time_range = 2000:6000;         % X-axis range for plotting
% ----------------------------------------------------------------


%% 1) Define a helper function to compute mean & SEM (excluding zeros)
computeMeanSEM = @(data) deal( ...
    nanmean(data(data~=0)), ...                       % mean of nonzero
    nanstd(data(data~=0)) / sqrt(sum(data~=0)) ...    % SEM of nonzero
);

%% 2) Compute mean & SEM for each timepoint in apr_x
num_timepoints = size(apr_x,2);
mean_apr_x = nan(1, num_timepoints);
sem_apr_x  = nan(1, num_timepoints);
for t = 1:num_timepoints
    [mean_apr_x(t), sem_apr_x(t)] = computeMeanSEM(apr_x(:,t));
end

%% 3) Repeat for apr_y, av_x, av_y
mean_apr_y = nan(1, num_timepoints);
sem_apr_y  = nan(1, num_timepoints);
for t = 1:num_timepoints
    [mean_apr_y(t), sem_apr_y(t)] = computeMeanSEM(apr_y(:,t));
end

mean_av_x = nan(1, num_timepoints);
sem_av_x  = nan(1, num_timepoints);
for t = 1:num_timepoints
    [mean_av_x(t), sem_av_x(t)] = computeMeanSEM(av_x(:,t));
end

mean_av_y = nan(1, num_timepoints);
sem_av_y  = nan(1, num_timepoints);
for t = 1:num_timepoints
    [mean_av_y(t), sem_av_y(t)] = computeMeanSEM(av_y(:,t));
end

%% 4) (Optional) Smooth the mean & SEM
sm_mean_apr_x = smooth(mean_apr_x, smoothing);
sm_sem_apr_x  = smooth(sem_apr_x,  smoothing);

sm_mean_apr_y = smooth(mean_apr_y, smoothing);
sm_sem_apr_y  = smooth(sem_apr_y,  smoothing);

sm_mean_av_x = smooth(mean_av_x, smoothing);
sm_sem_av_x  = smooth(sem_av_x,  smoothing);

sm_mean_av_y = smooth(mean_av_y, smoothing);
sm_sem_av_y  = smooth(sem_av_y,  smoothing);

%% 5) Plot PSTHs with shading for SEM
figure; 

% Define x-values shifted by 6000 so that x=6000 becomes time=0
xvals = (1:num_timepoints) - 6000;

% ---------------- Subplot 1: X-variables ------------------------
subplot(1,2,1); 
hold on;

% "apr_x": fill area
apr_x_upper = sm_mean_apr_x + sm_sem_apr_x*2;
apr_x_lower = sm_mean_apr_x - sm_sem_apr_x*2;
fill([xvals, fliplr(xvals)], ...
     [apr_x_upper', fliplr(apr_x_lower')], ...
     'b', 'FaceAlpha',0.3, 'EdgeColor','none');
% Plot the smoothed mean
plot(xvals, sm_mean_apr_x, 'b', 'LineWidth', 1.5);

% "av_x": fill area
av_x_upper = sm_mean_av_x + sm_sem_av_x*2;
av_x_lower = sm_mean_av_x - sm_sem_av_x*2;
fill([xvals, fliplr(xvals)], ...
     [av_x_upper', fliplr(av_x_lower')], ...
     'r', 'FaceAlpha',0.3, 'EdgeColor','none');
% Plot the smoothed mean
plot(xvals, sm_mean_av_x, 'r', 'LineWidth', 1.5);

% Shifted x-limits
xlim([time_range(1), time_range(end)] - 6000);
ylim([-0.5 4.5]);

xlim([-1300 0])
title('PSTH for X variables');
xlabel('Time (aligned so that 6000 = 0)');
ylabel('Mean value (\pm 2 SEM)');
box off;                 % remove box
set(gca, 'TickDir','out');   % ticks outward

% ---------------- Subplot 2: Y-variables ------------------------
subplot(1,2,2); 
hold on;

% "apr_y"
apr_y_upper = sm_mean_apr_y + sm_sem_apr_y*2;
apr_y_lower = sm_mean_apr_y - sm_sem_apr_y*2;
fill([xvals, fliplr(xvals)], ...
     [apr_y_upper', fliplr(apr_y_lower')], ...
     'b', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(xvals, sm_mean_apr_y, 'b-', 'LineWidth', 1.5);

% "av_y"
av_y_upper = sm_mean_av_y + sm_sem_av_y*2;
av_y_lower = sm_mean_av_y - sm_sem_av_y*2;
fill([xvals, fliplr(xvals)], ...
     [av_y_upper', fliplr(av_y_lower')], ...
     'r', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(xvals, sm_mean_av_y, 'r-', 'LineWidth', 1.5);

xlim([time_range(1), time_range(end)] - 6000);
ylim([-0.5 4.5]);
title('PSTH for Y variables');
xlabel('Time (aligned so that 6000 = 0)');
ylabel('Mean value (\pm 2 SEM)');
box off;
set(gca, 'TickDir','out');

xlim([-1300 0])
