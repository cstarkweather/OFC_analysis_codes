%% plot nullclines for 3 trial types
alphareward= 0.08;
alphapunish = 0.13;
noise = 80;
w_i = 0.08;
alpha = 22;
lambda = 2.5;
offset = -22;
valence= 0.5;
combos = [0 4; 7 4; 7 0]; % each pair is reward/punishment

figure

for i=1:3
    subplot(1,3,i)
mfsim_fitting_plot_nullclines_figure(combos(i,1),combos(i,2), 0.001, 1,alphareward,alphapunish,noise,lambda,w_i,alpha,valence,offset,0)
end

exportgraphics(gcf, '/results/figure4a.png');
disp('Saved /results/figure4a.png');


load('/data/posterior_params.mat')
% posteriors for approach and avoid state (these will depend on the noise level) are pre-computed in 'calculate_posteriors' and saved as 'posterior_params.mat'


%%high conflict example trajectory


figure
hold on
% for a particular reward/punishment contingency (here, reward = 7; punishment = 4), compute nullclines and the trajectory with some random gaussian noise on a single example trial
[trajx trajy]=mfsim_fitting_plot_nullclines_figure(7, 4, 0.001, 1,alphareward,alphapunish,noise,lambda,w_i,alpha,valence,offset,1);
hold on

smoothing_kernel = 5;
p = zeros(1,length(trajx)/smoothing_kernel);
epochs = smoothing_kernel: smoothing_kernel:length(trajx);

for i=1:length(epochs)-1
    p(i) = calculate_posterior(mean(trajx(epochs(i):epochs(i)+smoothing_kernel)),mean(trajy(epochs(i):epochs(i)+smoothing_kernel)),mu_x_appr,sd_x_appr,mu_y_appr,sd_y_appr,mu_x_av,sd_x_av,mu_y_av,sd_y_av,pi_appr,pi_av);
end

% once the posterior probability of a particular state stabilizes to favor one state over the other (p(state) = 0.9) for a certain amount of time (duration used in manuscript = 400ms), make a decision. make note of that decision time to just plot the example trajectory until that time
[decision_time, decision, num_state]=detectStatePreferences(p,400,0.9)

%% another way to plot data above over time (different timepoints shown with color gradient)
xdata = smooth(trajx(1:decision_time) , 20);
ydata = smooth(trajy(1:decision_time) , 20);

xdata = resample(xdata,800,1000);
ydata = resample(ydata,800,1000);

N = length(xdata);  % number of points in the trajectory
alphaVal = 0.2;

for i = 1:N

    thisColor = [1-i/N, i/N, 1];

scatter( xdata(i), ydata(i), ...
         10, ...                      % Marker size
         'o', ...                    % Marker shape
         'LineWidth', 2, ...         % Thicker outline
         'MarkerEdgeColor', thisColor, ...
         'MarkerFaceColor', thisColor, ...
         'MarkerEdgeAlpha', alphaVal, ...
         'MarkerFaceAlpha', alphaVal );
    hold on
end

axis equal;  % if appropriate
box off;
set(gca, 'TickDir', 'out');
xlim([-1.5 5.5])
ylim([-1.5 5.5])
title('high conflict example trajectory')

exportgraphics(gcf, '/results/figure4b.png');
disp('Saved /results/figure4b.png');





%%low conflict example trajectory

figure
hold on
% for a particular reward/punishment contingency (here, reward = 7; punishment = 0), compute nullclines and the trajectory with some random gaussian noise on a single example trial
[trajx trajy]=mfsim_fitting_plot_nullclines_figure(7, 0, 0.001, 1,alphareward,alphapunish,noise,lambda,w_i,alpha,valence,offset,1);
hold on

smoothing_kernel = 5;
p = zeros(1,length(trajx)/smoothing_kernel);
epochs = smoothing_kernel: smoothing_kernel:length(trajx);

for i=1:length(epochs)-1
    p(i) = calculate_posterior(mean(trajx(epochs(i):epochs(i)+smoothing_kernel)),mean(trajy(epochs(i):epochs(i)+smoothing_kernel)),mu_x_appr,sd_x_appr,mu_y_appr,sd_y_appr,mu_x_av,sd_x_av,mu_y_av,sd_y_av,pi_appr,pi_av);
end

% once the posterior probability of a particular state stabilizes to favor one state over the other (p(state) = 0.9) for a certain amount of time (duration used in manuscript = 400ms), make a decision. make note of that decision time to just plot the example trajectory until that time
[decision_time, decision, num_state]=detectStatePreferences(p,400,0.9)

%% another way to plot data above over time (different timepoints shown with color gradient)
xdata = smooth(trajx(1:decision_time) , 20);
ydata = smooth(trajy(1:decision_time) , 20);

xdata = resample(xdata,800,1000);
ydata = resample(ydata,800,1000);

N = length(xdata);  % number of points in the trajectory
alphaVal = 0.2;

for i = 1:N

    thisColor = [1-i/N, i/N, 1];

scatter( xdata(i), ydata(i), ...
         10, ...                      % Marker size
         'o', ...                    % Marker shape
         'LineWidth', 2, ...         % Thicker outline
         'MarkerEdgeColor', thisColor, ...
         'MarkerFaceColor', thisColor, ...
         'MarkerEdgeAlpha', alphaVal, ...
         'MarkerFaceAlpha', alphaVal );
    hold on
end

axis equal;  % if appropriate
box off;
set(gca, 'TickDir', 'out');
xlim([-1.5 5.5])
ylim([-1.5 5.5])
title('low conflict example trajectory')

exportgraphics(gcf, '/results/figure4c.png');
disp('Saved /results/figure4c.png');




%% now compute p_approach for all trial types


%reward and punishment contingencies
r = [2 2 2 4 4 4 7 7 7 0 0];
p = [0 2 4 0 2 4 0 2 4 2 4];

n_trials = 20;
decision = zeros(1,n_trials);
p_approach = zeros(1,length(r));
conflict = zeros(1,length(r));
avg_num_states = zeros(1,length(r));
num_states = zeros(length(r),n_trials);
trajectory = struct();


for trial_type = 1:11
    reward = r(trial_type);
    punishment = p(trial_type);
    for n = 1:n_trials

        [trajx,trajy] = mfsim_fitting_plot_nullclines_figure(reward,punishment,0, 0,alphareward,alphapunish,noise,lambda,w_i,alpha,valence,offset,1);

        trajectory(trial_type).trajectories_x(n,:) = trajx;
        trajectory(trial_type).trajectories_y(n,:) = trajy;

        posteriors = zeros(1,length(trajx));

        for i=1:length(trajx)
            posteriors(i) = calculate_posterior(trajx(i),trajy(i),mu_x_appr,sd_x_appr,mu_y_appr,sd_y_appr,mu_x_av,sd_x_av,mu_y_av,sd_y_av,pi_appr,pi_av);
        end

        trajectory(trial_type).posterior(n,:) = posteriors;

    end


end
%% Figure out decisions and num_states, plot psychometric curves
% Example data: row vector of posterior probabilities from 0.1

trial_type=struct();
time_stable = 400;
thresh = 0.9;

noise = randn(1,n_trials);

for i = 1:11
    num_states = zeros(1,n_trials);
    decisions = zeros(1,n_trials);
    decision_times = zeros(1,n_trials);
    for n = 1:n_trials
        [decision_time, decision, num_state] = detectStatePreferences(smooth(trajectory(i).posterior(n,:),30),time_stable,thresh);
        decision_times(n) = decision_time;
        decisions(n) = decision;
        num_states(n) = num_state;
    end
    trial_type(i).decisions = decisions;
    trial_type(i).decision_times = decision_times;
    trial_type(i).num_states = num_states;
    p_approach(i) = sum(decisions)/length(decisions);
    prob = p_approach(i);
    if prob ==0 || prob ==1
        conflict(i) = 0;
    else
        conflict(i) = -prob * log2(prob) - (1 - prob) * log2(1 - prob);
    end
end


% Assuming r, p, trial_type, trialID, and decision are defined
reward_sizes = [0 2 4 7];

r = [2 2 2 4 4 4 7 7 7 0 0];
p = [0 2 4 0 2 4 0 2 4 2 4];


%% relationship between conflict and num_states visited per trial
for i = 1:11
    avg_decision_time(i) = sum(trial_type(i).decision_times)/length(trial_type(i).decision_times);
    avg_num_states(i)= sum(trial_type(i).num_states)/length(trial_type(i).num_states);
    avg_transition_rate(i) = sum(trial_type(i).num_states)/sum(trial_type(i).decision_times)*1000;
end

%% num_states visited per trial versus conflict, and state transition rate versus conflict
figure
subplot(1,2,1)
% 1) Plot raw data as black dots
plot(conflict, avg_num_states, 'k.', 'MarkerSize', 20);
box off;
set(gca, 'TickDir', 'out');
hold on;

% 2) Compute correlation (Pearson by default)
[rho, pval] = corr(conflict(:), avg_num_states(:), 'Type','Pearson');

% 3) Fit a linear regression line (y = p(1)*x + p(2))
p = polyfit(conflict(:), avg_num_states(:), 1);

% 4) To get a neat line from left to right, sort the x-values
[xSorted, idx] = sort(conflict);
yFit = polyval(p, xSorted);

% 5) Plot best-fit line in black
plot(xSorted, yFit, 'k-', 'LineWidth', 2);

% 6) Display correlation results in the title (or use text(...) somewhere in the axes)
title(sprintf('r = %.3f, p = %.3g', rho, pval));
xlabel('conflict');
ylabel('avg num states');
xlim([0 1])


subplot(1,2,2)
plot(conflict, avg_transition_rate, 'k.', 'MarkerSize', 20);
box off;
set(gca, 'TickDir', 'out');
hold on;

% 2) Compute correlation (Pearson by default)
[rho, pval] = corr(conflict(:), avg_transition_rate(:), 'Type','Pearson');

% 3) Fit a linear regression line (y = p(1)*x + p(2))
p = polyfit(conflict(:), avg_transition_rate(:), 1);

% 4) To get a neat line from left to right, sort the x-values
[xSorted, idx] = sort(conflict);
yFit = polyval(p, xSorted);

% 5) Plot best-fit line in black
plot(xSorted, yFit, 'k-', 'LineWidth', 2);

% 6) Display correlation results in the title (or use text(...) somewhere in the axes)
title(sprintf('r = %.3f, p = %.3g', rho, pval));
xlabel('conflict');
ylabel('transition rate');
xlim([0 1])

exportgraphics(gcf, '/results/figure4d.png');
disp('Saved /results/figure4d.png');