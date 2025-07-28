%% use decoding approach to identify posteriors on each trial
%run this entire code first, which performs regression for the decoder and saves posteriors
%then run Fig_5
%this code presents stepwise AIC, which produces similar results to Lasso
clear sliding_epochs;
spacing = 100;start_sliding = 4000;length_epoch = 300;
last_sliding = 5000 - length_epoch;num_epochs = (last_sliding - start_sliding)/spacing + 1;
sliding_epochs= zeros(num_epochs,2);
sliding_epochs(1,:) = [start_sliding (start_sliding +length_epoch-1)];

for i=2:num_epochs
    sliding_epochs(i,:) = [(sliding_epochs(i-1,1)+spacing) (sliding_epochs(i-1,2)+ spacing)];
end

%% find the best decoding electrodes and epoch for each subject

% Build a minimum AIC decoding model for each epoch, and obtain accuracy using bootstrapping
num_subjects = 6;
subject_array = 1:6;
num_bootstraps = 500;

% Initialize storage for results
all_subjects_accuracies = cell(num_subjects, num_epochs);
all_subjects_mean_accuracies = zeros(num_subjects, num_epochs);
all_subjects_shuffled_accuracies = cell(num_subjects, num_epochs);
all_subjects_mean_shuffled_accuracies = zeros(num_subjects, num_epochs);
all_subjects_selected_predictors = cell(num_subjects, num_epochs);

train_ratio = 0.9;

% Loop through subjects
for j = 1:num_subjects
    n = subject_array(j);
    fprintf('Processing Subject %d...\n', n);

    % Filter indices based on decisions
    filter = find(subject(n).decision<2);
    num_electrodes = length(subject(n).electrode);

    % Loop through epochs
    for i = 1:num_epochs
        % Initialize the predictor matrix for this epoch
        valid_trials = 1:length(filter);
        predictors = zeros(length(valid_trials), num_electrodes);

                for e = 1:num_electrodes
                    if subject(n).electrode(e).olf_coor>-1
                    signal = sum(subject(n).electrode(e).trigger(2).high_gamma_mat(filter(valid_trials), sliding_epochs(i, 1):sliding_epochs(i, 2)), 2) / length_epoch;
                    predictors(:, e) = signal;
                    else
                    end
                end

                % Response variable
                response_variable = subject(n).decision(filter(valid_trials));

        % % Perform forward selection using AIC
        stepwise_model = stepwiseglm(predictors, response_variable, ...
            'constant', ...  % Start with no predictors
            'upper', 'linear', ...  % Allow all predictors in the final model
            'Distribution', 'binomial', ...
            'Criterion', 'AIC', ...
            'Verbose', 0); % Set to 1 for detailed output

        % Store the indices of predictors included in the model for this epoch
        selected_predictor_indices = stepwise_model.Formula.InModel;
        all_subjects_selected_predictors{j, i} = find(selected_predictor_indices);

        % Bootstrap to assess model accuracy
        bootstrap_accuracies = zeros(num_bootstraps, 1);
        bootstrap_shuffled_accuracies = zeros(num_bootstraps, 1);

        for b = 1:num_bootstraps
            % Split the data into training and testing sets
            cv = cvpartition(length(response_variable), 'Holdout', 1 - train_ratio);
            train_idx = training(cv);
            test_idx = test(cv);

            % Train the model on the training set with selected predictors
            trained_model = fitglm(predictors(train_idx, all_subjects_selected_predictors{n, i}), response_variable(train_idx), ...
                                   'linear', ...
                                   'Distribution', 'binomial');

            % Test the model on the testing set
            predicted_response = predict(trained_model, predictors(test_idx, all_subjects_selected_predictors{n, i}));
            binary_predicted_response = predicted_response > 0.5;

            % Calculate accuracy for this bootstrap
            bootstrap_accuracies(b) = sum(binary_predicted_response == response_variable(test_idx)) / length(response_variable(test_idx));

            % Shuffle the response variable and calculate shuffled accuracy
            shuffled_response = response_variable(randperm(length(response_variable)));
            shuffled_model = fitglm(predictors(train_idx, all_subjects_selected_predictors{n, i}), shuffled_response(train_idx), ...
                                    'linear', ...
                                    'Distribution', 'binomial');
            shuffled_predicted_response = predict(shuffled_model, predictors(test_idx, all_subjects_selected_predictors{n, i}));
            binary_shuffled_response = shuffled_predicted_response > 0.5;
            bootstrap_shuffled_accuracies(b) = sum(binary_shuffled_response == shuffled_response(test_idx)) / length(shuffled_response(test_idx));
        end

        % Store bootstrap accuracies and mean accuracy for this epoch
        all_subjects_accuracies{j, i} = bootstrap_accuracies;
        all_subjects_mean_accuracies(j, i) = mean(bootstrap_accuracies);
        all_subjects_shuffled_accuracies{j, i} = bootstrap_shuffled_accuracies;
        all_subjects_mean_shuffled_accuracies(j, i) = mean(bootstrap_shuffled_accuracies);
    end
end



% Calculate confidence intervals and plot
% Extract the time points (you can adjust the offset as needed)
timepoints = sliding_epochs(:,1)'-5000;
num_epochs = size(sliding_epochs,1);
% Initialize storage for standard errors
all_subjects_se_lower = zeros(num_subjects, num_epochs);
all_subjects_se_upper = zeros(num_subjects, num_epochs);
all_subjects_shuffled_se_lower = zeros(num_subjects, num_epochs);
all_subjects_shuffled_se_upper = zeros(num_subjects, num_epochs);

% Calculate standard errors for each subject and epoch
for n = 1:num_subjects
    for i = 1:num_epochs
        actual_accuracies = all_subjects_accuracies{n, i};
        shuffled_accuracies = all_subjects_shuffled_accuracies{n, i};

        % Calculate standard error for actual accuracies
        se_actual = std(actual_accuracies) / sqrt(num_bootstraps);
        all_subjects_se_lower(n, i) = mean(actual_accuracies) - se_actual*2;
        all_subjects_se_upper(n, i) = mean(actual_accuracies) + se_actual*2;

        % Calculate standard error for shuffled accuracies
        se_shuffled = std(shuffled_accuracies) / sqrt(num_bootstraps);
        all_subjects_shuffled_se_lower(n, i) = mean(shuffled_accuracies) - se_shuffled*2;
        all_subjects_shuffled_se_upper(n, i) = mean(shuffled_accuracies) + se_shuffled*2;
    end
end

% Plotting for all subjects
figure;
for subject_index = 1:num_subjects
    % Extract relevant data for plotting
    mean_actual_accuracy = all_subjects_mean_accuracies(subject_index, :);
    mean_shuffled_accuracy = all_subjects_mean_shuffled_accuracies(subject_index, :);
    actual_se_intervals = [all_subjects_se_lower(subject_index, :)', all_subjects_se_upper(subject_index, :)'];
    shuffled_se_intervals = [all_subjects_shuffled_se_lower(subject_index, :)', all_subjects_shuffled_se_upper(subject_index, :)'];

    % Create a subplot for the current subject
    subplot(2, 3, subject_index);
    hold on;

    % Plot actual accuracy with standard errors
    plot(timepoints, mean_actual_accuracy, 'b', 'LineWidth', 2);
    fill([timepoints, fliplr(timepoints)], ...
        [actual_se_intervals(:, 1)', fliplr(actual_se_intervals(:, 2)')], ...
        'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Plot shuffled accuracy with standard errors
    plot(timepoints, mean_shuffled_accuracy, 'r', 'LineWidth', 2);
    fill([timepoints, fliplr(timepoints)], ...
        [shuffled_se_intervals(:, 1)', fliplr(shuffled_se_intervals(:, 2)')], ...
        'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Formatting
    box off;
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'Color', 'w');
    xlabel('Time from decision (ms)');
    ylabel('Accuracy');
    title(['Subject ' num2str(subject_index)]);
    xlim([min(timepoints), max(timepoints)]);
    ylim([0.4, 1]);
    grid on;
    hold off;
    %xlim([-1000 -500])
end

% Add a main title for the entire figure
sgtitle('Decoding Accuracies Over Time with Standard Errors for All Subjects', 'FontSize', 16);

%
% Calculate mean and standard error for actual and shuffled accuracies across subjects
mean_actual = mean(all_subjects_mean_accuracies, 1);  % Mean across subjects (rows)
mean_shuffled = mean(all_subjects_mean_shuffled_accuracies, 1);  % Mean across subjects (rows)

std_actual = std(all_subjects_mean_accuracies, 0, 1);  % Standard deviation across subjects
std_shuffled = std(all_subjects_mean_shuffled_accuracies, 0, 1);  % Standard deviation across subjects

% Standard error: std / sqrt(number of subjects)
num_subjects = size(all_subjects_mean_accuracies, 1);
se_actual = std_actual / sqrt(num_subjects);  % Standard error for actual accuracy
se_shuffled = std_shuffled / sqrt(num_subjects);  % Standard error for shuffled accuracy

% Define time points (x-axis), assuming sliding_epochs(:,2) gives time points
time_points = sliding_epochs(:, 2) - 5000;  % Adjust time points as needed

% Create the figure
figure;
hold on;

% Plot mean actual accuracy with error bars
errorbar(time_points, mean_actual, se_actual, '-o', 'DisplayName', 'Actual Accuracy', 'LineWidth', 1.5, 'Color', 'b');

% Plot mean shuffled accuracy with error bars
errorbar(time_points, mean_shuffled, se_shuffled, '-o', 'DisplayName', 'Shuffled Accuracy', 'LineWidth', 1.5, 'Color', 'r');

% Add labels, title, and legend
xlabel('time from decision (ms)');
ylabel('Accuracy');
title('Mean Accuracy with Error Bars (Actual vs. Shuffled)');
legend('show');
box off
hold off;
set(gca, 'TickDir', 'out');
ylim([0.5 0.7])
% Define tick marks, ensuring that 0 is included
y_limits = ylim;
y_ticks = linspace(y_limits(1), y_limits(2), 3);  % Create 5 evenly spaced ticks
% Set the y-axis ticks
set(gca, 'YTick', y_ticks);
xlim([-550 50])

% save accuracies

save('/data/accuracies.mat','all_subjects_selected_predictors','all_subjects_accuracies',...
    'all_subjects_mean_accuracies','all_subjects_shuffled_accuracies','all_subjects_mean_shuffled_accuracies','sliding_epochs');


%%

clear trial_epochs;
spacing = 5;start_sliding = 1;length_trial_epoch = 10;
last_sliding = 5001 - length_trial_epoch;num_epochs1 = (last_sliding - start_sliding)/spacing + 1;
trial_epochs= zeros(num_epochs1,2);
trial_epochs(1,:) = [start_sliding (start_sliding +length_trial_epoch-1)];
for i=2:num_epochs1
    trial_epochs(i,:) = [(trial_epochs(i-1,1)+spacing) (trial_epochs(i-1,2)+ spacing) ];
end

% Go through every trial and compute posterior probabilities (unshuffled and shuffled)

num_shuffles = 100; % Number of shuffles for null distribution
for n=1:6

    best_epoch = find(all_subjects_mean_accuracies(n,:) == max(all_subjects_mean_accuracies(n,:)));

    % Get the predictors and response variable for the best epoch
    filter = find(subject(n).decision<2);

    % Response variable
    response_variable = subject(n).decision(filter);

    electrode_indices = all_subjects_selected_predictors{n,best_epoch};

    predictors = zeros(length(filter),length(electrode_indices));

    for i=1:length(electrode_indices)
        predictors(:,i) = sum(subject(n).electrode(electrode_indices(i)).trigger(2).high_gamma_mat(filter, sliding_epochs(best_epoch, 1):sliding_epochs(best_epoch, 2)), 2) / (sliding_epochs(best_epoch, 2)-sliding_epochs(best_epoch, 1)+1);
    end

    % Train the model using the reduced predictors
    glm_model = fitglm(predictors, response_variable, 'Distribution', 'binomial');

    % Initialize posterior probabilities
    num_trials = length(filter);
    num_epochs1 = size(trial_epochs, 1);
    posterior_probabilities = NaN(num_trials, num_epochs1);

    posterior_probabilities_shuffled = cell(1, num_shuffles);

    for j=1:num_shuffles
        posterior_probabilities_shuffled{j} = NaN(num_trials, num_epochs1);

        % Train the SHUFFLED model using the reduced predictors
        shuffled_response = response_variable(randperm(length(response_variable)));
        glm_model_shuffled = fitglm(predictors, shuffled_response, 'Distribution', 'binomial');

        % Predict posterior probabilities for actual and shuffled data
        for x = 1:num_trials
            for i = 1:num_epochs1
                if (5000 - trial_epochs(i, 1)) > subject(n).rt(filter(x))
                else

                    % Initialize the predictor matrix
                    new_predictors = zeros(1, length(electrode_indices));

                    % Sum neural data for each electrode
                    for b = 1:length(electrode_indices)
                        signal = sum(subject(n).electrode(electrode_indices(b)).trigger(2).high_gamma_mat(filter(x), trial_epochs(i, 1):trial_epochs(i, 2)), 2) / length_trial_epoch;
                        new_predictors(b) = signal;
                    end


                    if j==1
                        % Predict the posterior probabilities based on real data
                        posterior_probs = predict(glm_model, new_predictors);
                        posterior_probabilities(x, i) = posterior_probs; % Since we want P(decision == 1)

                        % Predict the posterior probabilities based on shuffled data
                        posterior_probs_shuffled = predict(glm_model_shuffled, new_predictors);
                        posterior_probabilities_shuffled{j}(x, i) = posterior_probs_shuffled; % Since we want P(decision == 1)
                    else
                        % Predict the posterior probabilities based on shuffled data
                        posterior_probs_shuffled = predict(glm_model_shuffled, new_predictors);
                        posterior_probabilities_shuffled{j}(x, i) = posterior_probs_shuffled; % Since we want P(decision == 1)
                    end

                end
            end
        end
    end
    savename = strcat(['/data/posteriors_' num2str(n)]);
    save(savename,'posterior_probabilities','posterior_probabilities_shuffled','glm_model')

    savename = strcat(['/data/chosen_indices_' num2str(n)]);
    save(savename,'electrode_indices')
end

%% Extract information about states
num_subjects = 6;
combined_conf_interval = cell(1, num_subjects);

for n=1:num_subjects
    loadname = strcat(['posteriors_' num2str(n) '.mat']);
    load(loadname)

    %calculate the combined confidence interval
    % Assume posterior_probabilities_shuffled is a cell array where each cell contains a matrix
    num_shuffles = length(posterior_probabilities_shuffled);

    % Initialize an empty array to store all data points
    all_data = [];

    % Loop through each shuffled posterior probabilities matrix
    for i = 1:num_shuffles
        data = posterior_probabilities_shuffled{i};
        % Flatten the matrix to a vector and concatenate with all_data
        all_data = [all_data; data(:)];
    end

    % Compute the 95% confidence interval for all combined data points
    combined_conf_interval{n} = prctile(all_data, [2.5 97.5]);

    filter = find(subject(n).decision<2);
    % Initialize arrays to store entry times, departure times, state identities, trial identities, state durations, and average posterior probabilities
    entry_times = [];
    departure_times = [];
    state_identity = [];
    trial_identity = [];
    reward = [];
    conflict = [];
    value = [];
    punishment = [];
    state_duration = [];
    transfer_time = [];
    average_posterior_prob = [];  % New array to store the average posterior probability
    entry_indices = [];
    departure_indices = [];
    decision = [];rt=[];

    % Loop over each trial
    for trial = 1:size(posterior_probabilities, 1)
        % Extract the posterior probabilities for the current trial
        posterior_probs = posterior_probabilities(trial, :);

        % Initialize a state variable
        states = zeros(size(posterior_probs));

        % Identify 1 states (posterior probability > upper bound of CI)
        states(posterior_probs > combined_conf_interval{n}(2)) = 1;

        % Identify -1 states (posterior probability < lower bound of CI)
        states(posterior_probs < combined_conf_interval{n}(1)) = -1;

        defined_times = find(~isnan(posterior_probs));

        % Handle the first state if it's already outside of CI
        if states(defined_times(1)) ~= 0
            entry_times = [entry_times; NaN];
            entry_indices = [entry_indices;NaN];
            state_identity = [state_identity; states(defined_times(1))];
            trial_identity = [trial_identity; filter(trial)];
            reward = [reward; subject(n).reward_trial(filter(trial))];
            punishment = [punishment; subject(n).punish_trial(filter(trial))];
            conflict = [conflict; subject(n).conflict_trial(filter(trial))];
            value = [value; subject(n).value_trial(filter(trial))];
            decision = [decision; subject(n).decision(filter(trial))];
            rt = [rt; subject(n).rt(filter(trial))];

            % Find the departure time
            departure_idx = find(states(defined_times(1)+1:end) ~= states(defined_times(1)), 1) + defined_times(1);
            if isempty(departure_idx)
                departure_times = [departure_times; NaN];
                departure_indices = [departure_indices; NaN];
                state_duration = [state_duration; NaN];
                transfer_time = [transfer_time; NaN];
                average_posterior_prob = [average_posterior_prob; NaN];
            else
                departure_times = [departure_times; trial_epochs(departure_idx, 1) - 5000];
                departure_indices = [departure_indices; departure_idx];
                state_duration = [state_duration; NaN];
                transfer_time = [transfer_time; NaN];
                average_posterior_prob = [average_posterior_prob; mean(posterior_probs(defined_times(1):departure_idx))];
            end
        end

        % Loop through states to identify entries and departures
        for i = defined_times(2:end)
            % Entry into state 1 or -1
            if states(i) ~= 0 && states(i-1) ~= states(i)
                % Record entry time
                entry_times = [entry_times; trial_epochs(i, 1) - 5000];
                entry_indices = [entry_indices;i];

                % Record the state identity
                state_identity = [state_identity; states(i)];

                % Record the trial identity
                trial_identity = [trial_identity; filter(trial)];
                reward = [reward; subject(n).reward_trial(filter(trial))];
                punishment = [punishment; subject(n).punish_trial(filter(trial))];
                conflict = [conflict; subject(n).conflict_trial(filter(trial))];
                value = [value; subject(n).value_trial(filter(trial))];
                decision = [decision; subject(n).decision(filter(trial))];
                rt = [rt; subject(n).rt(filter(trial))];

                % Find the departure time
                departure_idx = find(states(i+1:end) ~= states(i), 1) + i;
                if isempty(departure_idx)
                    % No departure found, assign NaN
                    departure_times = [departure_times; NaN];
                    departure_indices = [departure_indices; NaN];
                    state_duration = [state_duration; - entry_times(end)];
                    transfer_time = [transfer_time; NaN];
                    average_posterior_prob = [average_posterior_prob; mean(posterior_probs(i:end))];
                else
                    % Record departure time
                    departure_times = [departure_times; trial_epochs(departure_idx, 1) - 5000];
                    departure_indices = [departure_indices; departure_idx];
                    state_duration = [state_duration; departure_times(end) - entry_times(end)];
                    transfer_time = [transfer_time; NaN]; % Initialize transfer time, to be calculated next loop
                    average_posterior_prob = [average_posterior_prob; mean(posterior_probs(i:departure_idx))];
                end
            end

            % Calculate transfer time between states
            if length(departure_times) > 1 && length(entry_times) == length(departure_times)
                transfer_time(end) = entry_times(end) - departure_times(end-1);
            end
        end
    end

    % Combine the vectors into a single matrix if needed
    state_times_matrix = [entry_times, departure_times, state_identity, trial_identity, state_duration, transfer_time, average_posterior_prob];


    % store state transition times
    % Initialize vectors for storing transition times and types
    state_transition_times = [];
    transitiondeparture_times_store = [];
    transitionentry_times_store = [];
    transition_types = [];
    transition_trial_identity=[];
    reward_t=[];punishment_t=[];rt_t=[];
    conflict_t=[];value_t = [];decision_t = [];

    % Loop over each trial identity
    for trial = unique(trial_identity)'
        % Get the relevant indices for the current trial
        trial_indices = find(trial_identity == trial);

        % Extract relevant data for the current trial
        current_posterior_probs = posterior_probabilities(filter==trial, :);
        current_entry_times = entry_indices(trial_indices);
        current_departure_times = departure_indices(trial_indices);
        current_state_identity = state_identity(trial_indices);

        % Loop through the states in the current trial
        for i = 1:length(current_departure_times)

            crossing_idx=[];


            if ~isnan(current_departure_times(i)) && length(current_entry_times)>=(i+1)
                % Get the current state and its entry/departure times
                next_state_entry = current_entry_times(i+1);
                state_departure = current_departure_times(i);
                state_id = current_state_identity(i);
                next_state_id = current_state_identity(i+1);

                if current_posterior_probs(state_departure)>0.5
                    % Identify the index where the state crosses 0.5
                    crossing_idx = find(current_posterior_probs(state_departure:next_state_entry) < combined_conf_interval{n}(2), 1);
                else
                    % Identify the index where the state crosses 0.5
                    crossing_idx = find(current_posterior_probs(state_departure:next_state_entry) > combined_conf_interval{n}(1), 1);
                end

            end

            if ~isempty(crossing_idx)
                % Store the crossing time
                transition_time = state_departure + crossing_idx - 1;
                state_transition_times = [state_transition_times; trial_epochs(transition_time,1)-5000];

                % Check the state before and after the crossing
                if state_id == -1 && next_state_id == 1
                    transition_type = 1;  % 1--> -1
                elseif state_id == 1 && next_state_id == -1
                    transition_type = -1;  % -1 --> 1
                elseif state_id == 1 && next_state_id == 1
                    transition_type = 2;  % 1 --> 1
                elseif state_id == -1 && next_state_id == -1
                    transition_type = -2;  % -1 --> -1
                else
                    transition_type = NaN;  % No valid transition detected
                end

                % Store the transition type, entry, and departure times
                transition_types = [transition_types; transition_type];
                transitiondeparture_times_store = [transitiondeparture_times_store; trial_epochs(state_departure,1)-5000];
                transitionentry_times_store = [transitionentry_times_store; trial_epochs(next_state_entry,1)-5000];
                transition_trial_identity = [transition_trial_identity;trial];
                reward_t = [reward_t; subject(n).reward_trial(trial)];
                punishment_t = [punishment_t; subject(n).punish_trial(trial)];
                conflict_t = [conflict_t; subject(n).conflict_trial(trial)];
                value_t = [value_t; subject(n).value_trial(trial)];
                decision_t = [decision_t; subject(n).decision(trial)];
                rt_t = [rt_t; subject(n).rt(trial)];
            end
        end

    end

    % find full states. also have a flag for if a full state corresponds to the first state entered in a trial.
    % Initialize arrays to store the flags and halfway times
    full_state_flag = zeros(size(entry_times));
    first_full_state_flag = zeros(size(entry_times));
    last_state_flag = zeros(size(entry_times));  % Initialize last state flag
    halfway_time = nan(size(entry_times));

    % Loop over each state to determine if it's a full state, calculate halfway time, and identify last state
    for i = 1:length(entry_times)
        % Check if both entry and departure times are defined
        if ~isnan(entry_times(i)) && ~isnan(departure_times(i))
            full_state_flag(i) = 1;  % Mark as full state

            % Calculate the halfway time between entry and departure
            halfway_time(i) = entry_times(i) + (departure_times(i) - entry_times(i)) / 2;
        end

        % Check if this full state is the first in its trial
        if full_state_flag(i) == 1 && (i == 1 || trial_identity(i) ~= trial_identity(i-1))
            first_full_state_flag(i) = 1;
        end

        % Check if this full state is the last in its trial
        if full_state_flag(i) == 1 && (i == length(entry_times) || trial_identity(i) ~= trial_identity(i+1))
            last_state_flag(i) = 1;
        end
    end

    savename = strcat(['/data/state_metadata_' num2str(n) '.mat']);
    save(savename,'full_state_flag','halfway_time','first_full_state_flag','last_state_flag','transition_trial_identity',...
        'transitionentry_times_store','transitiondeparture_times_store','transition_types', 'rt','rt_t',...
        'entry_times','departure_times','state_identity','trial_identity','state_duration', 'transfer_time', 'average_posterior_prob',...
        'reward','punishment','reward_t','punishment_t','state_transition_times','conflict','conflict_t','value','value_t','decision','decision_t','combined_conf_interval');

end
