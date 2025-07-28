%% Calculate posterior probability null distributions
% determine state of simulated neural trajectories
%% experiment
alphareward= 0.08;
alphapunish = 0.13;
noise = 80;
w_i = 0.08;
alpha = 22;
lambda = 2.5;
offset = -22;
valence= 0.5;

clear decisions

%reward and punishment contingencies
r = [2 2 2 4 4 4 7 7 7 0 0];
p = [0 2 4 0 2 4 0 2 4 2 4];

n_trials = 100;


trial_type = 11; %trial type where there is no conflict and the data must generated around 'avoid' state
reward = r(trial_type);
punishment = p(trial_type);

trajectories_x_av = zeros(n_trials,5000);
trajectories_y_av = zeros(n_trials,5000);

for n = 1:n_trials

    [trajx,trajy] = mfsim_fitting_plot_nullclines_figure(reward,punishment,0, 0,alphareward,alphapunish,noise,lambda,w_i,alpha,valence,offset,1);

    trajectories_x_av(n,:) = trajx(1001:end);
    trajectories_y_av(n,:) = trajy(1001:end);


end


trial_type = 7; %trial type where there is no conflict and the data must generated around 'approach' state
reward = r(trial_type);
punishment = p(trial_type);


trajectories_x_appr = zeros(n_trials,5000);
trajectories_y_appr = zeros(n_trials,5000);

for n = 1:n_trials

    [trajx,trajy] = mfsim_fitting_plot_nullclines_figure(reward,punishment,0, 0,alphareward,alphapunish,noise,lambda,w_i,alpha,valence,offset,1);

    trajectories_x_appr(n,:) = trajx(1001:end);
    trajectories_y_appr(n,:) = trajy(1001:end);


end

% Flatten each matrix into one column vector
x_appr = trajectories_x_appr(:);
y_appr = trajectories_y_appr(:);
x_av   = trajectories_x_av(:);
y_av   = trajectories_y_av(:);

% Fit univariate Gaussians
mu_x_appr = mean(x_appr);  sd_x_appr = std(x_appr);
mu_y_appr = mean(y_appr);  sd_y_appr = std(y_appr);

mu_x_av   = mean(x_av);    sd_x_av   = std(x_av);
mu_y_av   = mean(y_av);    sd_y_av   = std(y_av);

% Priors
N_appr = numel(x_appr);
N_av   = numel(x_av);
pi_appr = N_appr / (N_appr + N_av);
pi_av   = 1 - pi_appr;   % or use 0.5, etc.


smoothed_trajx = smooth(trajx,1);
smoothed_trajy = smooth(trajy,1);

posterior_appr = zeros(1,length(smoothed_trajx));
posterior_av = zeros(1,length(smoothed_trajx));

for i=1:length(smoothed_trajx)
x_current = smoothed_trajx(i);y_current=smoothed_trajy(end);
obs = [x_current, y_current];

% Evaluate likelihoods
px_curr_appr = normpdf(x_current, mu_x_appr, sd_x_appr);
py_curr_appr = normpdf(y_current, mu_y_appr, sd_y_appr);
pXY_appr = px_curr_appr * py_curr_appr;  % naive independence

px_curr_av   = normpdf(x_current, mu_x_av,   sd_x_av);
py_curr_av   = normpdf(y_current, mu_y_av,   sd_y_av);
pXY_av   = px_curr_av * py_curr_av;

% Posterior probability state=1
posterior_appr(i) = (pXY_appr * pi_appr) / ...
    (pXY_appr * pi_appr + pXY_av * pi_av);
posterior_av(i) = 1 - posterior_appr(i);

end
save('posterior_params.mat','mu_x_appr','sd_x_appr','mu_y_appr','sd_y_appr','mu_x_av','sd_x_av','mu_y_av','sd_y_av','pi_appr','pi_av')
%% nice plot to show the distributions

subplot(2,2,1)
histogram(trajectories_x_appr)
xlim([-2 5])
box off;
set(gca, 'TickDir', 'out');
hold off


subplot(2,2,3)
histogram(trajectories_x_av)
xlim([-2 5])
box off;
set(gca, 'TickDir', 'out');
hold off

subplot(2,2,2)
histogram(trajectories_y_appr)
xlim([-2 5])
box off;
set(gca, 'TickDir', 'out');
hold off

subplot(2,2,4)
histogram(trajectories_y_av)
xlim([-2 5])
box off;
set(gca, 'TickDir', 'out');
hold off
%%
