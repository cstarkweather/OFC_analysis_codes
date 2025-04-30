% Define the range of subjects to loop over
subjectRange = 5; %one in Figure 1 is Subject 5

% Initialize an empty matrix to store results
% Columns: [SubjectID, ModelType (1 for Conflict, 2 for Reaction Time), R_squared, Beta, P_value]
results1 = [];results2 = [];

alpha_reward_store = [];
alpha_punish_store = [];
lambda_store = [];
mu_store = [];

% Loop through each subject
for n = subjectRange

    % Call the plotSubjectData function and capture the R^2, beta weights, and p-values
    [rho_conflict, pval_conflict, rho_rt, pval_rt,alpha_reward,alpha_punish,lambda,mu] = ...
        plotSubjectData(n, subject(n).reward_coeff, subject(n).punish_coeff, ...
                        subject(n).reward_trial, subject(n).punish_trial, ...
                        subject(n).decision, subject(n).conflict_trial_type, ...
                        subject(n).rt/1000, subject(n).reward_trial_type, subject(n).punishment_trial_type,subject(n).p_approach_trial_type,subject(n).conflict_trial');
    
    % Append the results for conflict vs. decision function (subplot 3)
    results1 = [results1; n, rho_conflict, pval_conflict];

    % Append the results for reaction time vs. decision function (subplot 4)
    results2 = [results2; n, rho_rt, pval_rt];
    alpha_reward_store = [alpha_reward_store alpha_reward];
    alpha_punish_store = [alpha_punish_store alpha_punish];
    mu_store = [mu_store mu];
    lambda_store=[lambda_store lambda];
    
end

% Display the results
disp('Results Matrix (SubjectID, ModelType, R_squared, Beta, P_value):');
disp(results1);
disp(results2);
%%
figure; hold on
jitter = 0.15;
for n = 1:6
    scatter(1- jitter/2 + rand * jitter,alpha_reward_store(n), 30, ...
        [n*0.15 n*0.15 n*0.15], 'filled');  % Set RGB color and use 'filled'
    hold on
    scatter(2- jitter/2 + rand * jitter,alpha_punish_store(n), 30, ...
        [n*0.15 n*0.15 n*0.15], 'filled');  % Set RGB color and use 'filled'
    hold on
    scatter(3- jitter/2 + rand * jitter,lambda_store(n), 30, ...
        [n*0.15 n*0.15 n*0.15], 'filled');  % Set RGB color and use 'filled'
    hold on
end
set(gca, 'Box', 'off', 'TickDir', 'out');
set(gcf, 'Color', 'w');
%%
alpha_reward_store
alpha_punish_store
lambda_store

%%
psychtemp = [0 0 1 1 0 1];
anova1(alpha_reward_store',psychtemp')
anova1(alpha_punish_store',psychtemp')
anova1(lambda_store',psychtemp')
%%
versiontemp = [0 1 1 0 0 0];
anova1(alpha_reward_store',versiontemp')
anova1(alpha_punish_store',versiontemp')
anova1(lambda_store',versiontemp')
%%
load('behavior.mat')
% Loop through each subject
for n = subjectRange

    % Call the plotSubjectData function and capture the R^2, beta weights, and p-values
    [rho_conflict, pval_conflict, rho_rt, pval_rt] = ...
        plotSubjectData(1, 1, 1, ...
                        reward_trial, punish_trial, ...
                        decision, conflict_trial_type, ...
                        reaction_time/1000, reward_trial_type, punishment_trial_type,p_approach_trial_type,conflict_trial');
    
    % Append the results for conflict vs. decision function (subplot 3)
    results1 = [results1; n, rho_conflict, pval_conflict];

    % Append the results for reaction time vs. decision function (subplot 4)
    results2 = [results2; n, rho_rt, pval_rt];
end

% Display the results
disp('Results Matrix (SubjectID, ModelType, R_squared, Beta, P_value):');
disp(results1);
disp(results2);

%%
for i=1:6
if results1(i,3)<0.05
    plot(1 + rand *0.1,results1(i,2),'k.')
    hold on
else
    plot(1 + rand *0.1,results1(i,2),'r.')
    hold on
end

if results2(i,3)<0.05
    plot(1.5 + rand *0.1,results2(i,2),'k.')
else
    plot(1.5 + rand *0.1,results2(i,2),'r.')
end
end

xlim([0.8 1.8])
ylim([- 1 1])
box off; set(gca, 'TickDir', 'out');
%% Adapted plotSubjectData function using prospect theory model
function [rho_conflict, pval_conflict, rho_rt, pval_rt,alpha_reward_opt,alpha_punish_opt,lambda_opt,mu_opt] = ...
    plotSubjectData(subjectID, reward_coeff, punish_coeff, rewards_trials, punishments_trials, ...
                    decision, conflict_all, reaction_time_all, reward_trial_types, punishment_trial_types, p_approach,conflict_trial)
    filter = find(decision<2);
    % Use the passed-in trial data for plotting
    rewards = rewards_trials(filter);
    punishments = punishments_trials(filter);
    decisions = decision(filter);
    %conflict = conflict_all(filter);
    reaction_times = reaction_time_all(filter);
    conflict_trial = conflict_trial(filter);
    
    % For grid-based estimation, use the unique reward and punishment sizes:
    reward_sizes = 0:1:7;
    punishment_sizes = 0:1:4;
    
    % Build a matrix (grid) for the observed p_approach values.
    % Rows correspond to punishment_sizes and columns to reward_sizes.
    mat = NaN(length(punishment_sizes), length(reward_sizes));
    for i = 1:length(punishment_sizes)
        for j = 1:length(reward_sizes)
            idx = find(punishment_trial_types == punishment_sizes(i) & reward_trial_types == reward_sizes(j));
            if ~isempty(idx)
                mat(i, j) = p_approach(idx);  % If more than one index exists, this uses the available data.
            end
        end
    end

    % Flatten the grid and remove NaNs (to build model fitting data)
    [X_grid, Y_grid] = meshgrid(reward_sizes, punishment_sizes);
    X_flat = X_grid(:);
    Y_flat = Y_grid(:);
    P_flat = mat(:);
    valid_idx = ~isnan(P_flat);
    X_flat = X_flat(valid_idx);
    Y_flat = Y_flat(valid_idx);
    P_flat = P_flat(valid_idx);

    %%% Prospect Theory Model Fitting Section
    % Known outcome probabilities (as provided):
    p_reward = 0.42;      % (0.6 * 0.7)
    p_punishment = 0.4;
    
    % Estimate parameters for the model using the grid data.
    % The cost function (costfun) computes the negative log-likelihood:
    %   params = [mu, alpha_reward, alpha_punish, lambda]
    objfun = @(params) costfun(params, X_flat, Y_flat, P_flat, p_reward, p_punishment);
    init_params = [1, 1, 1, 1];  % initial guesses for [mu, α_reward, α_punish, λ]
    options = optimset('Display','off','TolFun',1e-8);
    [opt_params, ~] = fminsearch(objfun, init_params, options);
    
    % Unpack optimized parameters:
    mu_opt = opt_params(1);
    alpha_reward_opt = opt_params(2);
    alpha_punish_opt = opt_params(3);
    lambda_opt = opt_params(4);
    
    % Compute the decision value (V) on the grid:
    V_grid = p_reward * (X_grid.^alpha_reward_opt) - p_punishment * lambda_opt * (Y_grid.^alpha_punish_opt);
    
    %%% Plotting Section
    % SUBPLOT 2: Plot the decision function (image) and overlay the decision boundary (V = 0)
    subplot(1, 4, 2)
    imagesc(reward_sizes, punishment_sizes, V_grid)
    colormap bone
    colorbar
    hold on;
    % The decision boundary is where V = 0 (p_pred = 0.5)
    contour(reward_sizes, punishment_sizes, V_grid, [0 0], 'LineColor', 'k', 'LineWidth', 2);
    set(gca, 'YDir', 'normal');
    xlabel('Reward size');
    ylabel('Punishment size');
    title('Decision boundary (V=0)');
    box off; set(gca, 'TickDir', 'out');
    hold off;

    % SUBPLOT 1: Scatter plot of individual trials (with jitter)
    subplot(1, 4, 1)
    rand_factor = 0.3;
    scatter(rewards(decisions==1) + rand_factor*randn(sum(decisions==1),1), ...
            punishments(decisions==1) + rand_factor*randn(sum(decisions==1),1), ...
            50, 'o', 'MarkerEdgeColor', [0.0909, 0.9091, 0.1500]);
    hold on;
    scatter(rewards(decisions==0) + rand_factor*randn(sum(decisions==0),1), ...
            punishments(decisions==0) + rand_factor*randn(sum(decisions==0),1), ...
            50, '+', 'MarkerEdgeColor', [1.0000, 0, 0.1500]);
    % Overlay the decision boundary from the grid (for reference)
    contour(reward_sizes, punishment_sizes, V_grid, [0 0], 'LineColor', 'k', 'LineWidth', 2);
    xlim([-0.5 7.5]); ylim([-0.5 4.5]);
    xlabel('Reward size'); ylabel('Punishment size');
    title('Trial data & decision boundary');
    box off; set(gca, 'TickDir', 'out');
    hold off;
    
    % SUBPLOT 3: Conflict vs. |Value|
    subplot(1, 4, 3)
    % Compute trial-level decision value using the prospect theory model
    % (Using individual trial offers, with the same estimated exponents)
    Value = p_reward * (reward_trial_types.^alpha_reward_opt) - p_punishment * lambda_opt * (punishment_trial_types.^alpha_punish_opt);
    %decision_function_conflict = abs(Val_apr_conflict);
    % Plot the data points
    plot(abs(Value), conflict_all, 'k.', 'MarkerSize', 20);
    hold on;
    % Fit a linear regression model and plot the best-fit line
    [rho_conflict, pval_conflict] = corr(abs(Value)', conflict_all');
    mdl_conflict = fitlm(abs(Value), conflict_all);
    beta_conflict = mdl_conflict.Coefficients.Estimate(2);
    x_fit = linspace(min(abs(Value)), max(abs(Value)), 100);
    y_fit = mdl_conflict.Coefficients.Estimate(1) + beta_conflict*x_fit;
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    text(max(x_fit), max(y_fit), sprintf('p = %.3f, rho = %.3f', pval_conflict, rho_conflict), ...
         'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    xlabel('|Decision Value|'); ylabel('Conflict');
    title('Conflict vs. |Value|');
    box off; set(gca, 'TickDir', 'out');
    hold off;

    % SUBPLOT 4: Reaction Time vs. |Value|
    subplot(1, 4, 4)
    % Compute trial-level decision value for reaction time analysis
    Val_apr_rt = p_reward * (rewards.^alpha_reward_opt) - p_punishment * lambda_opt * (punishments.^alpha_punish_opt);
    decision_function_rt = abs(Val_apr_rt);
    scatter(decision_function_rt + 0.01*randn(length(decision_function_rt),1), reaction_times, 5, 'k.');
    hold on;
    [rho_rt, pval_rt] = corr(decision_function_rt(decisions < 2), reaction_times(decisions < 2));
    mdl_rt = fitlm(decision_function_rt, reaction_times);
    beta_rt = mdl_rt.Coefficients.Estimate(2);
    x_fit_rt = linspace(min(decision_function_rt), max(decision_function_rt), 100);
    y_fit_rt = mdl_rt.Coefficients.Estimate(1) + beta_rt*x_fit_rt;
    plot(x_fit_rt, y_fit_rt, 'r-', 'LineWidth', 2);
    text(max(x_fit_rt), max(y_fit_rt), sprintf('p = %.3f, rho = %.3f', pval_rt, rho_rt), ...
         'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    xlabel('|Decision Value|'); ylabel('Reaction time (s)');
    title('RT vs. |Value|');
    box off; set(gca, 'TickDir', 'out');
    hold off;

end

%% Cost Function for Prospect Theory Model
function err = costfun(params, X, Y, Pobs, p_reward, p_punishment)
    % params = [mu, alpha_reward, alpha_punish, lambda]
    mu = params(1);
    alpha_reward = params(2);
    alpha_punish = params(3);
    lambda = params(4);
    
    % Compute the approach value for each observation:
    V = p_reward * (X.^alpha_reward) - p_punishment * lambda * (Y.^alpha_punish);
    
    % Predicted probability via logistic function:
    p_pred = 1./(1+exp(-mu*V));
    
    % Clip probabilities to avoid log(0):
    epsilon = 1e-10;
    p_pred = min(max(p_pred, epsilon), 1-epsilon);
    
    % Negative log-likelihood:
    LL = sum(Pobs.*log(p_pred) + (1-Pobs).*log(1-p_pred));
    err = -LL;  % Minimization target
end
