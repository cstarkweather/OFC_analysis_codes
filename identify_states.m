%% use decoding approach to identify posteriors on each trial
%run this entire code first, which saves posteriors
%then run Fig_5
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



%% Calculate confidence intervals and plot
% Extract the time points (you can adjust the offset as needed)
timepoints = sliding_epochs(:,1)'-5000;  % Adjust this based on your time reference
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

%%
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

%% save accuracies

save('accuracies.mat','all_subjects_selected_predictors','all_subjects_accuracies',...
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
%% Go through every trial and compute posterior probabilities (unshuffled and shuffled)

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
        predictors(:,i) = sum(subject(n).electrode(electrode_indices(i)).trigger(3).high_gamma_mat(filter, sliding_epochs(best_epoch, 1):sliding_epochs(best_epoch, 2)), 2) / (sliding_epochs(best_epoch, 2)-sliding_epochs(best_epoch, 1)+1);
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
                        signal = sum(subject(n).electrode(electrode_indices(b)).trigger(3).high_gamma_mat(filter(x), trial_epochs(i, 1):trial_epochs(i, 2)), 2) / length_trial_epoch;
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
    savename = strcat(['posteriors_' num2str(n)]);
    save(savename,'posterior_probabilities','posterior_probabilities_shuffled','glm_model')

    savename = strcat(['chosen_indices_' num2str(n)]);
    save(savename,'electrode_indices')
end
%%
% Define a colormap for the subjects
subject_colors = [
    205,  92,  150;  % Dusty Pink 
    144, 238,  180;  % Light Neutral Green
    147, 112,  219;  % Milkshake Purple
    130, 100,  100;  % Taupe
    255, 153,  102;  % Creamsicle Orange
    255, 253,  209; % Cream White
] / 255;  % Normalize to [0, 1]

% Generate the figure
figure(1);
standard_font_size = 20; 
numbers_font_size = 20;

% Plot the ROC values
subplot('Position', [0.05, 0.1, 0.45, 0.8]); % Adjusted width to make space for other plots
imagesc(sorted_roc_tc);  % Ensure this data has appropriate values
hold on;
colormap pink;
title('ROC values', 'FontSize', standard_font_size);
xlabel('Time (bins)', 'FontSize', standard_font_size);
ylabel('Electrodes', 'FontSize', standard_font_size);
set(gca, 'FontSize', numbers_font_size);  % Set font size for axis numbers
set(gca, 'TickDir', 'out');  % Set tick marks to face out

% Check for clusters and add cluster lines
if ~isnan(pc_roc)
    clustlines = nan(3, nclust);
    for i = 1:nclust
        temp = sum(clustlabel <= i) + 0.5;
        clustlines(1:2, i) = [temp; temp];
    end
    clustlines = clustlines(:);
    xlines = cat(1, repmat(xlim', [1, nclust]), nan(1, nclust));
    xlines = xlines(:);
    plot(xlines, clustlines, 'r');
end

% Initialize the color block for electrode inclusion
electrode_inclusion_color_block = ones(length(sorted_electrode_ID), 1, 3);  % Initialize as white [1, 1, 1]
chosen_indices=cell(1,5);
% Loop through each subject and determine electrode inclusion
for n = 1:num_subjects
    indices_name = strcat(['chosen_indices_' num2str(n) '.mat']);
    load(indices_name);
    chosen_indices{n} = electrode_indices;

    % Find the indices of electrodes that belong to the current subject
    subject_indices = find(sorted_subject_ID == n);
    
    % Loop through the subject's electrodes and mark included electrodes as black
    for idx = 1:length(subject_indices)
        electrode_id = sorted_electrode_ID(subject_indices(idx));
        
        if ismember(electrode_id, chosen_indices{n})
            % Set color to black [0, 0, 0] if included
            electrode_inclusion_color_block(subject_indices(idx), 1, :) = [0, 0, 0];
        end
    end
end

% Plot the electrode inclusion color block
subplot('Position', [0.51, 0.1, 0.08, 0.8]);  % Adjusted position and width
imagesc(electrode_inclusion_color_block);
title('Decoding', 'FontSize', standard_font_size);
yticks([]);  % Remove Y-axis tick marks
xticks([]);  % Remove X-axis tick marks
set(gca, 'FontSize', numbers_font_size);  % Set font size for axis numbers
set(gca, 'TickDir', 'out');  % Set tick marks to face out
colormap(gca, [1 1 1; 0 0 0]);  % Set colormap to black and white

% Create a color block for the coordinates using the jet colormap
jet_colormap = jet(256);  % Jet colormap

% Normalize the coordinates to the range [1, 256] for indexing into the colormap
norm_coordinates = round(rescale(sorted_coordinates, 1, 256));

% Create the color block for the coordinates
coordinate_color_block = jet_colormap(norm_coordinates, :);

% Plot the coordinate color block
subplot('Position', [0.61, 0.1, 0.08, 0.8]);  % Adjusted position and width
imagesc(reshape(coordinate_color_block, length(sorted_coordinates), 1, 3));
title('Coordinates', 'FontSize', standard_font_size);
yticks([]);  % Remove Y-axis tick marks
xticks([]);  % Remove X-axis tick marks
set(gca, 'TickDir', 'out');  % Set tick marks to face out

% Add coordinate labels directly onto the color block
for i = 1:5:length(sorted_coordinates)
    text(0.5, i, num2str(sorted_coordinates(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numbers_font_size, 'Color', 'k');
end

% Create a color block for the subjects
subject_color_block = zeros(length(sorted_subjects), 1, 3);  % Initialize the color block
for i = 1:length(folders)
    subject_indices = find(strcmp(sorted_subjects, folders{i}));
    for j = 1:length(subject_indices)
        subject_color_block(subject_indices(j), 1, :) = subject_colors(i, :);
    end
end

% Plot the subject color block
subplot('Position', [0.71, 0.1, 0.08, 0.8]);  % Adjusted position and width
imagesc(subject_color_block);
title('Subjects', 'FontSize', standard_font_size);
yticks([]);  % Remove Y-axis tick marks
xticks([]);  % Remove X-axis tick marks
set(gca, 'TickDir', 'out');  % Set tick marks to face out

% Ensure the entire figure fits into the view
set(gcf, 'Position', [100, 100, 1500, 600]);  % Adjust the figure size as needed

%%

n = 3;
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
combined_conf_interval = prctile(all_data, [2.5 97.5]);

% Display the confidence intervals
disp('95% Confidence Interval for All Shuffled Posterior Probabilities:');
disp(combined_conf_interval);

% proportion of values outside CI

% Flatten the matrix to a vector
temp = posterior_probabilities(~isnan(posterior_probabilities));
posterior_probabilities_vector = temp(:);

% Count the number of elements outside the confidence interval
num_outside = sum(posterior_probabilities_vector < combined_conf_interval(1) | posterior_probabilities_vector > combined_conf_interval(2));

% Calculate the total number of elements
total_elements = numel(posterior_probabilities_vector);

% Calculate the proportion of elements outside the confidence interval
proportion_outside = num_outside / total_elements
%% Histograms comparing unshuffled vs shuffled posterior probability distribution

binwidth = 0.02;
fontSize = 14;  % Define font size

% Plot for unshuffled posterior probabilities
subplot(2,1,1)
histogram(posterior_probabilities_vector, 'BinWidth', binwidth)
title('Posterior Probability Distribution - Unshuffled', 'FontSize', fontSize)
xlabel('Posterior Probability', 'FontSize', fontSize)
set(gca, 'TickDir', 'out', 'Color', 'w')  % Set tick direction and background color

% Plot for shuffled posterior probabilities
subplot(2,1,2)
histogram(all_data, 'BinWidth', binwidth)
title('Posterior Probability Distribution - Shuffled', 'FontSize', fontSize)
xlabel('Posterior Probability', 'FontSize', fontSize)
set(gca, 'TickDir', 'out', 'Color', 'w')  % Set tick direction and background color

