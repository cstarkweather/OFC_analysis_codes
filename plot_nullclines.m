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
for i=1:3
    subplot(1,3,i)
mfsim_fitting_plot_nullclines_figure(combos(i,1),combos(i,2), 0.001, 1,alphareward,alphapunish,noise,lambda,w_i,alpha,valence,offset,0)
end

%% shows how activity can bounce between approach(green) and avoid(red) states

mfsim_fitting_plot_nullclines_figure(4,2, 0.001, 1,alphareward,alphapunish,noise,lambda,w_i,alpha,valence,offset,0)

hold on
load('posterior_params.mat') % posteriors for approach and avoid state are computed in 'calculate_posteriors'
[trajx trajy]=mfsim_fitting_plot_nullclines_figure(4, 2, 0.001, 0,alphareward,alphapunish,noise,lambda,w_i,alpha,valence,offset,1);

smoothing_kernel = 5;
p = zeros(1,length(trajx)/smoothing_kernel);
epochs = smoothing_kernel: smoothing_kernel:length(trajx);

for i=1:length(epochs)-1
    p(i) = calculate_posterior(mean(trajx(epochs(i):epochs(i)+smoothing_kernel)),mean(trajy(epochs(i):epochs(i)+smoothing_kernel)),mu_x_appr,sd_x_appr,mu_y_appr,sd_y_appr,mu_x_av,sd_x_av,mu_y_av,sd_y_av,pi_appr,pi_av);
    plot(mean(trajx(epochs(i):epochs(i)+smoothing_kernel)),mean(trajy(epochs(i):epochs(i)+smoothing_kernel)),'.','Markersize',10,'Color',[1-p(i) p(i) p(i)*0.5 ])
    drawnow
    hold on
end


%% another way to plot data above over time (different timepoints shown with color gradient)
cut = 1500;
xdata = smooth(trajx(1:cut) , 20);
ydata = smooth(trajy(1:cut) , 20);

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

%% Smooth posteriors and determine state
posteriors = p;
plot(posteriors)
hold on
plot(smooth(posteriors,50),'Linewidth',2)
[decision_time, decision, num_state]=detectStatePreferences(posteriors,400,0.9)

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

figure;
hold on;
p_approach_avg = p_approach;
for i = 1:length(reward_sizes)
    reward_size = reward_sizes(i);
    x_data = p(r == reward_size);
    y_data = p_approach_avg(r == reward_size);

    % Define the logistic function
    logistic_function = @(x, m, b) 1 ./ (1 + exp(-(m * x + b)));

    % Define the objective function for fitting
    objective_function = @(params) sum((y_data - logistic_function(x_data, params(1), params(2))).^2);

    % Initial guess for parameters
    initial_params = [1, 1];

    % Find the parameters that minimize the objective function using fminsearch
    params = fminsearch(objective_function, initial_params);
    color = [1 - (i-1) / length(reward_sizes), (i-1) / length(reward_sizes), 1];

    % Plot the data

    scatter(x_data, y_data, 'DisplayName', 'Data', 'MarkerEdgeColor', color);

    % Plot the logistic curve
    x_range = linspace(min(x_data), max(x_data), 100);
    y_fit = logistic_function(x_range, params(1), params(2));
    plot(x_range, y_fit, 'DisplayName', 'Logistic Fit', 'Color', color);
end

ylim([0 1.2]);
box off;
set(gca, 'TickDir', 'out');
xlabel('# Bombs', 'FontSize', 20);
ylabel('% Approached', 'FontSize', 20);

% Turn off title
title('');

hold off;

for i = 1:11
    avg_decision_time(i) = sum(trial_type(i).decision_times)/length(trial_type(i).decision_times);
    avg_num_states(i)= sum(trial_type(i).num_states)/length(trial_type(i).num_states);
    avg_transition_rate(i) = sum(trial_type(i).num_states)/sum(trial_type(i).decision_times)*1000;
end
avg_decision_time
%% relationship between conflict and num_states visited per trial
for i = 1:11
    avg_decision_time(i) = sum(trial_type(i).decision_times)/length(trial_type(i).decision_times);
    avg_num_states(i)= sum(trial_type(i).num_states)/length(trial_type(i).num_states);
    avg_transition_rate(i) = sum(trial_type(i).num_states)/sum(trial_type(i).decision_times)*1000;
end

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
