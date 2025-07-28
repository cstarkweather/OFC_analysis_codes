%% get the 'behavior.mat' file by running 'extract_behavioral_variables' in the data analysis code
load('behavior.mat')
    % Call the plotSubjectData function and capture the R^2, beta weights, and p-values
    [rho_conflict, pval_conflict, rho_rt, pval_rt] = ...
        plotSubjectData(reward_coeff, punish_coeff, ...
                        reward_trial, punish_trial, ...
                        decision, conflict_trial_type, ...
                        reaction_time, reward_trial_type, punishment_trial_type,p_approach_trial_type);


%% stats
mean(results2(:,3))
%% Modified plotSubjectData function to return R^2, beta, and p-values
function [rho_conflict, pval_conflict, rho_rt, pval_rt] = ...
    plotSubjectData(reward_coeff, punish_coeff, rewards_trials, punishments_trials, ...
                    decision, conflict_all, reaction_time_all, reward_trial_types, punishment_trial_types,p_approach)
    
    % Depending on the condition flag, select the appropriate data
    rewards = rewards_trials;
    punishments = punishments_trials;
    decisions = decision;
    conflict = conflict_all;
    reaction_times = reaction_time_all;
    reward_trials = reward_trial_types;
    punishment_trials = punishment_trial_types;

    % Define the reward and punishment sizes
    reward_sizes = 0:1:7;
    punishment_sizes = 0:1:4;

    % Create a grid of reward and punishment values
    %[X, Y] = meshgrid(reward_sizes, punishment_sizes);


mat= NaN(5,8);

    % Fill the matrix with your data
    for i = 1:5 % cycle through different punishment sizes
        for j = 1:8
            index = find(punishment_trial_types == punishment_sizes(i) & reward_trial_types == reward_sizes(j));
            if ~isempty(index)
                mat(i, j) = p_approach(index);
            end
        end
    end

    % Flatten the matrices
    [X, Y] = meshgrid(reward_sizes, punishment_sizes);
    X_flat = X(:);
    Y_flat = Y(:);
    P_flat = mat(:);

    % Remove NaN values
    valid_idx = ~isnan(P_flat);
    X_flat = X_flat(valid_idx);
    Y_flat = Y_flat(valid_idx);
    P_flat = P_flat(valid_idx);


% Fit the logistic regression model
% Use 'binomial' as the distribution type
[b, dev_full, stats] = glmfit([X_flat Y_flat], P_flat, 'binomial', 'link', 'logit', 'constant', 'off');

% b(1) is reward_coeff and b(2) is punish_coeff


    % Compute the value of the decision function (Val_apr - Val_av)
    Val_apr = reward_coeff * X + punish_coeff * Y;
    decision_function = Val_apr;

    % Plot the decision function as an image (subplot 2)
    subplot(1, 4, 2)
    imagesc(reward_sizes, punishment_sizes, decision_function)
    colormap bone
    colorbar
    hold on

    % Plot the decision boundary where p_approach = 50% (U_apr = 0)
    contour(reward_sizes, punishment_sizes, decision_function, [0 0], 'LineColor', 'k', 'LineWidth', 2);

    % Ensure the y-axis goes from bottom (smallest) to top (largest)
    set(gca, 'YDir', 'normal');

    % Add labels and title
    xlabel('Reward size')
    ylabel('Punishment size')
    title('Decision boundary for p_{approach} = 50%')

    % Set plot properties
    box off;
    set(gca, 'TickDir', 'out');
    set(gcf, 'Color', 'w');
    hold off

    % Plot scatter plot for decision trials (subplot 1)
    subplot(1, 4, 1)
    rand_factor = 1;

    % Scatter plot for trials where the decision was 1
    scatter(rewards(decisions == 1) + rand_factor * rand(sum(decisions == 1), 1) - rand_factor / 2, ...
            punishments(decisions == 1) + rand_factor * rand(sum(decisions == 1), 1) - rand_factor / 2, ...
            50, 'o', 'MarkerEdgeColor', [0.0909, 0.9091, 0.1500]); % Marker size is set to 50
    hold on;

    % Scatter plot for trials where the decision was 0
    scatter(rewards(decisions == 0) + rand_factor * rand(sum(decisions == 0), 1) - rand_factor / 2, ...
            punishments(decisions == 0) + rand_factor * rand(sum(decisions == 0), 1) - rand_factor / 2, ...
            50, '+', 'MarkerEdgeColor', [1.0000, 0, 0.1500]); % Marker size is set to 50

    % Plot the decision boundary where p_approach = 50% (U_apr = 0)
    contour(reward_sizes, punishment_sizes, decision_function, [0 0], 'LineColor', 'k', 'LineWidth', 2);
    xlim([-0.5 7.5]);
    ylim([-0.5 4.5]);

    % Add labels and title
    xlabel('Reward size');
    ylabel('Punishment size');
    title('Decision boundary for p_{approach} = 50%');

    % Set plot properties
    box off;
    set(gca, 'TickDir', 'out');
    set(gcf, 'Color', 'w');
    hold off;

    % Plot conflict vs. decision function (subplot 3)
    subplot(1, 4, 3)
    % Compute the value of the decision function (Val_apr - Val_av)
    Val_apr = reward_coeff * reward_trials + punish_coeff * punishment_trials;
    decision_function_conflict = abs(Val_apr);

    % Plot the data with larger marker size
    plot(decision_function_conflict, conflict, 'k.', 'MarkerSize', 20);
    hold on;

    % Fit a linear model
    [rho_conflict pval_conflict] = corr(decision_function_conflict', conflict');
    mdl_conflict = fitlm(decision_function_conflict, conflict);
    beta_conflict = mdl_conflict.Coefficients.Estimate(2); % Slope (beta weight)

    % Plot the best fit line
    x_fit = linspace(min(decision_function_conflict), max(decision_function_conflict), 100);
    y_fit = mdl_conflict.Coefficients.Estimate(1) + beta_conflict * x_fit;
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    text(max(x_fit), max(y_fit), sprintf('p = %.3f, rho = %.3f', pval_conflict, rho_conflict), ...
        'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    hold off;

    % Add labels and title for subplot 3
    xlabel('|Value|');
    ylabel('Conflict');
    title('Conflict vs. |Value|');
    box off;
    set(gca, 'TickDir', 'out');
    set(gcf, 'Color', 'w');

    % Plot reaction time vs. decision function (subplot 4)
    subplot(1, 4, 4)
    scatter_factor = 0.5;
    Val_apr_rt = reward_coeff * rewards + punish_coeff * punishments;
    decision_function_rt = abs(Val_apr_rt);

    % Plot the data with larger marker size
    plot(decision_function_rt + rand(length(decision_function_rt),1) * scatter_factor, reaction_times, 'k.', 'MarkerSize', 5);
    hold on;
    [rho_rt pval_rt] = corr(decision_function_rt(decision <2), reaction_times(decision <2));
    % Fit a linear model for reaction time
    mdl_rt = fitlm(decision_function_rt, reaction_times);
    beta_rt = mdl_rt.Coefficients.Estimate(2); % Slope (beta weight)
    p_value_rt = mdl_rt.Coefficients.pValue(2); % p-value for the slope
    r_squared_rt = mdl_rt.Rsquared.Ordinary; % R-squared value

    % Plot the best fit line
    x_fit_rt = linspace(min(decision_function_rt), max(decision_function_rt), 100);
    y_fit_rt = mdl_rt.Coefficients.Estimate(1) + beta_rt * x_fit_rt;
    plot(x_fit_rt, y_fit_rt, 'r-', 'LineWidth', 2);
    text(max(x_fit_rt), max(y_fit_rt), sprintf('p = %.3f, rho = %.3f', pval_rt, rho_rt), ...
        'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    hold off;

    % Add labels and title for subplot 4
    xlabel('|Value|');
    ylabel('Reaction time (ms)');
    title('Reaction time vs. |Value|');
    box off;
    set(gca, 'TickDir', 'out');
    set(gcf, 'Color', 'w');
end
