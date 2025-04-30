%% Extract information about states
num_subjects = 6;
combined_conf_interval = cell(1, num_subjects);

for n=1:num_subjects
    loadname = strcat(['posteriors_100_' num2str(n) '.mat']);
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

    savename = strcat(['state_metadata_' num2str(n) '.mat']);
    save(savename,'full_state_flag','halfway_time','first_full_state_flag','last_state_flag','transition_trial_identity',...
        'transitionentry_times_store','transitiondeparture_times_store','transition_types', 'rt','rt_t',...
        'entry_times','departure_times','state_identity','trial_identity','state_duration', 'transfer_time', 'average_posterior_prob',...
        'reward','punishment','reward_t','punishment_t','state_transition_times','conflict','conflict_t','value','value_t','decision','decision_t','combined_conf_interval');

end

%% plot accuracies over time
clear error_real
clear error_shuffled
load('accuracies.mat')
n = 4;
num_bootstraps = 500;
for i=1:size(sliding_epochs,1)
    error_real(i) = std(all_subjects_accuracies{n,i})/sqrt(500);
    error_shuffled(i) = std(all_subjects_shuffled_accuracies{n,i})/sqrt(500);
    plot(sliding_epochs(i,1),all_subjects_mean_accuracies(n,i),'.','Markersize',20,'Color',[0 0 0])
    hold on
    plot(sliding_epochs(i,1),all_subjects_mean_shuffled_accuracies(n,i),'.','Markersize',20,'Color',[0.5 0.5 0.5])
    hold on
end

errorbar(sliding_epochs(:,1)',all_subjects_mean_accuracies(n,:),error_real,'Color',[0 0 0],'CapSize', 0)
hold on
errorbar(sliding_epochs(:,1)',all_subjects_mean_shuffled_accuracies(n,:),error_shuffled,'Color',[0.5 0.5 0.5],'CapSize', 0)
hold on

yline(sum(subject(n).decision==1)/(sum(subject(n).decision==1) + sum(subject(n).decision==0)),'r--')

set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'w');
set(gcf, 'Color', 'w');
xlim([4250 4700])
%
%% example single trials - posterior probabilities over time - Figure 3E-F
figure(2)
n=5;
loadname = strcat(['posteriors_' num2str(n) '.mat']);
load(loadname)

loadname = strcat(['state_metadata_' num2str(n) '.mat']);
load(loadname)

loadname = strcat(['chosen_indices_' num2str(n) '.mat']);
load(loadname)
load('trial_epochs.mat')

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

which_trial = 66;


filter = find(subject(n).decision<2);
defined_bins = ~isnan(posterior_probabilities(which_trial,:));

% Plot the posterior probability for approach
plot(trial_epochs(defined_bins)-5000, posterior_probabilities(which_trial, defined_bins), 'LineWidth', 2, 'DisplayName', 'Approach Probability');
hold on;

% Plot the posterior probability for avoid (1 - approach)
plot(trial_epochs(defined_bins)-5000, 1 - posterior_probabilities(which_trial, defined_bins), 'LineWidth', 2, 'DisplayName', 'Avoid Probability');
hold on;

% Plot confidence interval lines
yline(combined_conf_interval{n}(1), 'LineWidth', 2);
hold on;
yline(combined_conf_interval{n}(2), 'LineWidth', 2);

% Label the plot and adjust formatting
ylim([0 1]);
xlabel('Time from cue (ms)', 'FontSize', 12);
ylabel('Posterior Probability (Approach)', 'FontSize', 12);
%legend({'Approach Probability', 'Avoid Probability'}, 'FontSize', 12);
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'w');
set(gcf, 'Color', 'w');
title(['Trial: ' num2str(which_trial) ', Decision: ' num2str(subject(n).decision(filter(which_trial)))], 'FontSize', 14);

hold on
xline(entry_times(trial_identity==filter(which_trial) & state_identity>0),'r')
xline(departure_times(trial_identity==filter(which_trial)  & state_identity>0),'r--')
xline(entry_times(trial_identity==filter(which_trial) & state_identity<0),'b')
xline(departure_times(trial_identity==filter(which_trial)  & state_identity<0),'b--')



set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'w');
set(gcf, 'Color', 'w');
%% examining state occupancy over time
figure(2)
clear sliding_epochs;
%for the data aligned to decision
spacing = 5; start_sliding = -2000; length_epoch = 10; 
last_sliding = 0 - length_epoch; num_epochs = (last_sliding - start_sliding) / spacing + 1;


bins = zeros(num_epochs, 2);
bins(1, :) = [start_sliding (start_sliding + length_epoch - 1)];

for i = 2:num_epochs
    bins(i, :) = [(bins(i-1, 1) + spacing) (bins(i-1, 2) + spacing)];
end

for condition = 1
    % Reset accumulation variables for each condition
    total_p_state_1 = zeros(size(bins, 1), 1);
    total_p_state_2 = zeros(size(bins, 1), 1);
    total_valid_bins = zeros(size(bins, 1), 1);

    % Initialize arrays to accumulate p_state values for each subject
    p_state_1_subjects = zeros(size(bins, 1), 6);
    p_state_2_subjects = zeros(size(bins, 1), 6);

    for n = 1:6
        state_metadata_name = strcat(['state_metadata_' num2str(n) '.mat']);
        load(state_metadata_name)


        
        switch condition
            case 1
                filter = find(subject(n).decision == 1);
        end

        p_state = zeros(size(bins));
        valid_bins = zeros(size(bins, 1), 1);

        for i = 1:length(filter)
            entry_times_trial = entry_times(trial_identity == filter(i));
            departure_times_trial = departure_times(trial_identity == filter(i));
            state_identity_trial = state_identity(trial_identity == filter(i));

            for j = 1:size(bins, 1)
                bin_ind = find(entry_times_trial < bins(j, 1) & (departure_times_trial > bins(j, 2) | isnan(departure_times_trial)));
                
                % Only consider bins where the response time is greater than the bin end
                if subject(n).rt(filter(i)) > -(bins(j, 2))
                    valid_bins(j) = valid_bins(j) + 1;
                end
                
                for a = 1:length(bin_ind)
                    if state_identity_trial(bin_ind(a)) == 1
                        p_state(j, 1) = p_state(j, 1) + 1;
                    else
                        p_state(j, 2) = p_state(j, 2) + 1;
                    end
                end
            end
        end

        % Normalize p_state values by the number of valid bins
        p_state(:, 1) = p_state(:, 1) ./ valid_bins;
        p_state(:, 2) = p_state(:, 2) ./ valid_bins;

        smoothing_bin=10;
        time_vec1=bins(1:end , 1);
        %plot(time_vec1, smooth(p_state(:, 1),10),'Color', [0.15 0.8 / condition 0.2 / condition,0.2], 'LineWidth', 1)
        %hold on
        %plot(time_vec1, smooth(p_state(:, 2),10),'Color', [0.9 0.15 * condition 0.15 * condition,0.2], 'LineWidth', 1)
        %hold on

        % Store the p_state values for each subject
        p_state_1_subjects(:, n) = p_state(:, 1);
        p_state_2_subjects(:, n) = p_state(:, 2);

        % Accumulate p_state and valid_bins for average calculation across subjects
        total_p_state_1 = total_p_state_1 + p_state(:, 1) .* valid_bins;
        total_p_state_2 = total_p_state_2 + p_state(:, 2) .* valid_bins;
        total_valid_bins = total_valid_bins + valid_bins;
    end

    % Calculate averages of p_state across subjects
    avg_p_state_1 = total_p_state_1 ./ total_valid_bins;
    avg_p_state_2 = total_p_state_2 ./ total_valid_bins;

    % Calculate standard error across subjects (n = 6)
    se_p_state_1 = std(p_state_1_subjects, 0, 2) ./ sqrt(6);
    se_p_state_2 = std(p_state_2_subjects, 0, 2) ./ sqrt(6);

    % Smooth the averages and standard error
    smoothing_bin = 5;
    smoothed_avg_1 = smooth(avg_p_state_1, smoothing_bin);
    smoothed_avg_2 = smooth(avg_p_state_2, smoothing_bin);
    smoothed_se_1 = smooth(se_p_state_1, smoothing_bin);
    smoothed_se_2 = smooth(se_p_state_2, smoothing_bin);

    % Time vector for plotting
    time_vec = bins(1:end - smoothing_bin, 1);

    % Plot the average p_state values with bold lines
    plot(time_vec, smoothed_avg_1(1:end - smoothing_bin), 'Color', [0.15 0.8 / condition 0.2 / condition], 'LineWidth', 2);
    hold on;
    plot(time_vec, smoothed_avg_2(1:end - smoothing_bin), 'Color', [0.9 0.15 * condition 0.15 * condition], 'LineWidth', 2);

    % Plot the standard error as shaded areas
    fill([time_vec; flipud(time_vec)], ...
         [smoothed_avg_1(1:end - smoothing_bin) - smoothed_se_1(1:end - smoothing_bin); flipud(smoothed_avg_1(1:end - smoothing_bin) + smoothed_se_1(1:end - smoothing_bin))], ...
         [0.15 0.8 / condition 0.2 / condition], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([time_vec; flipud(time_vec)], ...
         [smoothed_avg_2(1:end - smoothing_bin) - smoothed_se_2(1:end - smoothing_bin); flipud(smoothed_avg_2(1:end - smoothing_bin) + smoothed_se_2(1:end - smoothing_bin))], ...
         [0.9 0.15 * condition 0.15 * condition], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % Adjust plot appearance
    set(gca, 'TickDir', 'out');
    box off;
    set(gcf, 'Color', 'w');
    xlabel('Time from trial onset (ms)', 'FontSize', 20);
    ylabel('p(state)', 'FontSize', 24);
    hold on;

    xlim([-1400 10]); %for data aligned to time of decision
    %xlim([0 1500]); %for data aligned to time of cue

end

ylim([0 0.3]);

%% Noisy ramp vs discrete attractor?
pos_coding_e=[];
neg_coding_e=[];
pos_coding_n=[];
neg_coding_n=[];
for n = 1:6
    posteriors_name = strcat(['posteriors_' num2str(n) '.mat']);
    load(posteriors_name,'glm_model')
    electrodes_name = strcat(['chosen_indices_' num2str(n) '.mat']);
    load(electrodes_name)
    electrode_indices = electrode_indices';
    beta_weights = glm_model.Coefficients{2:end, 1};
    pos_coding_e = [pos_coding_e;electrode_indices(beta_weights>0)];
    neg_coding_e = [neg_coding_e;electrode_indices(beta_weights<0)];
    pos_coding_n = [pos_coding_n;ones(sum(beta_weights>0),1)*n];
    neg_coding_n = [neg_coding_n;ones(sum(beta_weights<0),1)*n];
end

colorflag = 0;
signal_bounds = [-300 300];
baseline_bounds = [4800 5000];
window = signal_bounds(2) - signal_bounds(1) + 1;
duration_threshold = 100;

smoothing_kernel = 1;ID = 1;

x2 = round(signal_bounds(1):smoothing_kernel:signal_bounds(2));

figure;
hold on;
[B,I] = sort(plotorder);

n_set = pos_coding_n;
e_set = pos_coding_e;

for_psth = zeros(length(e_set),signal_bounds(2) - signal_bounds(1)+1);


for j = 1:6

    for a = 1:length(e_set)
        e = e_set(a);
        n = n_set(a);
        state_name = strcat(['state_metadata_' num2str(n) '.mat']);
        load(state_name)
        bound1 = -rt*0.33;
        bound2 = -rt*0.67;
        trialset = [...
             {find(state_identity==-1 & entry_times>bound1)},...
             {find(state_identity==-1 & entry_times<bound1 & entry_times>bound2)},...
             {find(state_identity==-1 & entry_times<bound2)},...
             {find(state_identity==1 & entry_times<bound2)},...
             {find(state_identity==1 & entry_times<bound1 & entry_times>bound2)},...
             {find(state_identity==1 & entry_times>bound1)},
             ];


        indices = trialset{j};
        temp=[];
        for i=1:length(indices)
            if (entry_times(indices(i))+5000>abs(signal_bounds(1)))
                baseline = [subject(n).electrode(e).trigger(1).high_gamma_mat(trial_identity(indices(i)), baseline_bounds(1):baseline_bounds(2))];
                avg_baseline = sum(baseline, 2) / size(baseline, 2);
                temp = [temp;(subject(n).electrode(e).trigger(t).high_gamma_mat(trial_identity(indices(i)), (5000+entry_times(indices(i))+signal_bounds(1)):(5000+entry_times(indices(i))+signal_bounds(2))) - avg_baseline) ./ avg_baseline];
            else
            end
        end

        for_psth(a,:)=sum(temp)/size(temp,1);

    end
    y = smooth(sum(for_psth) / size(for_psth, 1), smoothing_kernel);
    dy = zeros(length(y),1);
    for c = 1:length(x2) - 1
        dy((c - 1) * smoothing_kernel + 1:c * smoothing_kernel) = std(sum(for_psth(:, (c - 1) * smoothing_kernel + 1:c * smoothing_kernel),2)/smoothing_kernel)/sqrt(size(for_psth,1));
    end

    smooth_y = y(1:smoothing_kernel:end+1-smoothing_kernel)';
    smoothed_dy = smooth(dy,smoothing_kernel);
    if colorflag ==0
        fill([(signal_bounds(1):signal_bounds(2)), fliplr((signal_bounds(1):signal_bounds(2)))], [(y - smoothed_dy)'  fliplr((y + smoothed_dy)')],[1 - (j - 1) / (length(trialset) - 0.9), (j - 1) / (length(trialset) - 0.9), 0.15], 'linestyle', 'none','FaceAlpha', 0.6);
        hold on
        plot((signal_bounds(1):signal_bounds(2))', y, 'Linewidth',2,'Color',[0.7*(1 - (j - 1) / (length(trialset) - 1)) 0.5*(j - 1) / (length(trialset) - 1)  0 1])
        hold on
    else

        fill([(signal_bounds(1):signal_bounds(2)), fliplr((signal_bounds(1):signal_bounds(2)))], [(y - smoothed_dy)'  fliplr((y + smoothed_dy)')], [1 - (j-1) / (length(trialset)-1), (j -1)/ (length(trialset)-1), 1], 'linestyle', 'none','FaceAlpha', 0.6);
        hold on
        plot((signal_bounds(1):signal_bounds(2))', y, 'k','Linewidth',2)
        hold on
    end

        for_anova{j} = for_psth;
end


% Adjust plot appearance
set(gca, 'TickDir', 'out');
box off;
set(gcf, 'Color', 'w');
xlabel('Time from state switch (milliseconds)', 'FontSize', 14);
ylabel('High gamma amplitude', 'FontSize', 14);
xline(0)
xlim([-20 120])

%% anova FOR POSITIVE STATE SWITCHES
anova2([sum(for_anova{4}(:,300:375),2),sum(for_anova{5}(:,300:375),2),sum(for_anova{6}(:,300:375),2)])

%% &&&&&&&&&&&&&&&&&&&&&&&&&&&&& STATISTICS IN THE TEXT 

%% state durations
state_durations = [];
for n = 1:6
    load(strcat(['state_metadata_' num2str(n) '.mat']))
    subplot(6,1,n)
    histogram(state_duration,'binwidth',10)
    box off
    set(gca, 'TickDir', 'out', 'Color', 'w')  % Set tick direction and background color
    xlim([0 300])
    state_durations=[state_durations;state_duration];
end

histogram(state_durations,'binwidth',10)
box off
set(gca, 'TickDir', 'out', 'Color', 'w')  % Set tick direction and background color
xlim([0 400])
mean(state_durations(~isnan(state_durations)))
std(state_durations(~isnan(state_durations)))
%% number of electrodes chosen
for n= 1:6
        load(strcat(['chosen_indices_' num2str(n) '.mat']))
        num_chosen(n) = length(electrode_indices);
end

mean(num_chosen)
std(num_chosen)
%% proportion of posterior probability outside of the 95% CI
proportion_outside0 = zeros(1,num_subjects);
for n = 1:6
    load(strcat(['posteriors_' num2str(n) '.mat']))
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
    proportion_outside0(n) = num_outside / total_elements;
end

mean(proportion_outside0)
std(proportion_outside0)
%% confirm proportion outside not just due to decoding window

load('accuracies.mat')

proportion_outside = zeros(1,num_subjects);
for n = 1:6
    load(strcat(['posteriors_' num2str(n) '.mat']))
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

    % proportion of values outside CI

    % Flatten the matrix to a vector
    % determine optimal decoding window
    best_epoch = sliding_epochs(find(all_subjects_mean_accuracies(n,:) == max(all_subjects_mean_accuracies(n,:))),1);

    temp1 = posterior_probabilities(:,1:round(best_epoch/5));
    posterior_probabilities_vector = temp1(~isnan(temp1));

    % Count the number of elements outside the confidence interval
    num_outside = sum(posterior_probabilities_vector < combined_conf_interval(1) | posterior_probabilities_vector > combined_conf_interval(2));

    % Calculate the total number of elements
    total_elements = numel(posterior_probabilities_vector);

    % Calculate the proportion of elements outside the confidence interval
    proportion_outside(n) = num_outside / total_elements;
end


mean(proportion_outside)
std(proportion_outside)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUPPLEMENTARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 4
%% state durations histogram - supplementary figure 4b
state_durations = [];
for n = 1:6
    state_metadata_name = strcat(['state_metadata_' num2str(n) '.mat']);
    load(state_metadata_name)
    state_durations = [state_durations;state_duration(~isnan(state_duration))];
end
histogram(state_durations,'binwidth',10)
mean(state_durations)
std(state_durations)

set(gca, 'TickDir', 'out'); % Tick marks outside
box off; % Remove the box

xlim([0 400])


%% state transitions for high and low conflict trials  - Figure 4A
figure(1)

num_subjects = 6;
transition_rate = zeros(num_subjects,size(trialset,2));

for n = 1:6

    trialset = {find(subject(n).conflict_trial'>0.1 & subject(n).decision<2)...
        find(subject(n).conflict_trial'<0.1 & subject(n).decision<2)};

    state_metadata_name = strcat(['state_metadata_' num2str(n) '.mat']);
    load(state_metadata_name)

    for a=1:length(trialset)

        filter = trialset{a};
        total_transitions = 0;total_time = 0;

        for i=1:length(filter)
            transition_times_trial = transitionentry_times_store(transition_trial_identity==filter(i));
            transition_type_trial = transition_types(transition_trial_identity==filter(i));
            %total_transitions = total_transitions + length(transition_type_trial);
            total_transitions = total_transitions + sum(transition_type_trial>-3);
            total_time = total_time + subject(n).rt(filter(i));

        end

        transition_rate(n,a) = total_transitions/total_time * 1000;
    end

end

plot(transition_rate','k-','Linewidth',2)
%xlabel('low conflict         high conflict', 'FontSize', 16);
ylabel('state transitions/second', 'FontSize', 16);
set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'w');
set(gcf, 'Color', 'w');
xlim([0.5 2.5])


[h p ci stats] = ttest(transition_rate(:,1),transition_rate(:,2))
%%
%% state transitions for high and low conflict trials  - Figure 4A
figure(1)

num_subjects = 6;
transition_rate = zeros(num_subjects,size(trialset,2));
conflict_all = [];transition_rate_all=[];subject_ids = [];
for n = 1:6

    clear trialset;
    unique_trial_types = unique(subject(n).trial_type_trial);
    for i = 1:length(unique_trial_types)
        trialset{i} = find(subject(n).trial_type_trial==unique_trial_types(i) & subject(n).decision<2);
    end

    state_metadata_name = strcat(['state_metadata_' num2str(n) '.mat']);
    load(state_metadata_name)

    for a=1:length(trialset)

        filter = trialset{a};
        total_transitions = 0;total_time = 0;

        for i=1:length(filter)
            transition_times_trial = transitionentry_times_store(transition_trial_identity==filter(i));
            transition_type_trial = transition_types(transition_trial_identity==filter(i));
            %total_transitions = total_transitions + length(transition_type_trial);
            total_transitions = total_transitions + sum(transition_type_trial>-3);
            total_time = total_time + subject(n).rt(filter(i));

        end

        transition_rate(n,a) = total_transitions/total_time * 1000;
    end

plot(subject(n).conflict_trial_type,transition_rate(n,unique_trial_types),'.','Markersize',20)
hold on
[h_ind(n) p_ind(n)] = corr(subject(n).conflict_trial_type',transition_rate(n,unique_trial_types)');

subject_ids = [subject_ids ones(1,length(unique_trial_types))*n];
conflict_all = [conflict_all subject(n).conflict_trial_type];
transition_rate_all = [transition_rate_all transition_rate(n,unique_trial_types)];

end

h_ind
p_ind

%[h p] = corr(conflict_all',transition_rate_all')
%% Compute transition rate for all subjects
figure(1)

num_subjects = 6;
transition_rate = zeros(num_subjects, size(trialset,2));
conflict_all = []; transition_rate_all = []; subject_ids = [];

% Create a color gradient from black to light grey for each subject.
colors = zeros(num_subjects, 3);
for n = 1:num_subjects
    % Linear interpolation: subject 1 is black ([0,0,0]), subject 6 is light grey ([0.8,0.8,0.8])
    c = (n-1) * (0.8/(num_subjects-1));
    colors(n,:) = [c, c, c];
end

% Set the jitter magnitude (adjust as needed)
jitterMagnitude = 0.01;

for n = 1:num_subjects
    clear trialset;
    unique_trial_types = unique(subject(n).trial_type_trial);
    for i = 1:length(unique_trial_types)
        trialset{i} = find(subject(n).trial_type_trial == unique_trial_types(i) & subject(n).decision < 2);
    end

    state_metadata_name = strcat(['state_metadata_' num2str(n) '.mat']);
    load(state_metadata_name)

    for a = 1:length(trialset)
        filter = trialset{a};
        total_transitions = 0;
        total_time = 0;
        for i = 1:length(filter)
            transition_times_trial = transitionentry_times_store(transition_trial_identity == filter(i));
            transition_type_trial = transition_types(transition_trial_identity == filter(i));
            total_transitions = total_transitions + sum(transition_type_trial > -3);
            total_time = total_time + subject(n).rt(filter(i));
        end
        transition_rate(n, a) = total_transitions / total_time * 1000;
    end

    % Use the original x-values (conflict) and add jitter
    x_original = subject(n).conflict_trial_type;
    x_jittered = x_original + jitterMagnitude * randn(size(x_original));
    
    % Use scatter to plot the jittered points with 40% opacity.
    scatter(x_jittered, transition_rate(n, unique_trial_types), 36, colors(n,:), 'filled', 'MarkerFaceAlpha', 0.4)
    hold on
    [h_ind(n), p_ind(n)] = corr(x_original', transition_rate(n, unique_trial_types)');
    
    subject_ids = [subject_ids, ones(1, length(unique_trial_types)) * n];
    conflict_all = [conflict_all, x_original];
    transition_rate_all = [transition_rate_all, transition_rate(n, unique_trial_types)];
    
    % Compute the best-fit line using the original (non-jittered) data.
    p_coeff = polyfit(x_original, transition_rate(n, unique_trial_types), 1);
    x_fit = linspace(min(x_original), max(x_original), 100);
    y_fit = polyval(p_coeff, x_fit);
    
    % Plot the best-fit line with the same subject-specific color.
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', colors(n,:));
end

xlabel('Conflict')
ylabel('Transition Rate')
title('Best Fit Lines for Each Subject')

%% ask whether conflict modulates state transition rate
tbl = table(subject_ids', conflict_all', transition_rate_all',...
    'VariableNames', {'Subject', 'Conflict', 'TransitionRate'});

tbl.Subject = categorical(tbl.Subject);

% Define a mixed-effects model:
% The fixed effect is 'Conflict', and we include a random intercept for each subject.
lme = fitlme(tbl, 'TransitionRate ~ Conflict + Subject');

% View the model results:
disp(lme)


%% calculate number of states visited as a function of conflict

figure(1)

num_subjects = 6;
num_states = zeros(num_subjects,size(trialset,2));
conflict_all = [];total_states_all=[];subject_ids = [];
for n = 1:6

    clear trialset;
    unique_trial_types = unique(subject(n).trial_type_trial);
    for i = 1:length(unique_trial_types)
        trialset{i} = find(subject(n).trial_type_trial==unique_trial_types(i) & subject(n).decision<2);
    end

    state_metadata_name = strcat(['state_metadata_' num2str(n) '.mat']);
    load(state_metadata_name)

    for a=1:length(trialset)

        filter = trialset{a};
        total_states = 0;

        for i=1:length(filter)
            entry_times_trial = entry_times(trial_identity==filter(i));
            total_states = total_states + length(entry_times_trial);

        end

        num_states(n,a) = total_states/length(filter);
    end

plot(subject(n).conflict_trial_type,num_states(n,unique_trial_types),'.','Markersize',20)
hold on
[h_ind(n) p_ind(n)] = corr(subject(n).conflict_trial_type',num_states(n,unique_trial_types)');

subject_ids = [subject_ids ones(1,length(unique_trial_types))*n];
conflict_all = [conflict_all subject(n).conflict_trial_type];
total_states_all = [total_states_all num_states(n,unique_trial_types)];

end

h_ind
p_ind


%% nicer plot of all the slopes for all subjects (same data as above)
figure(1)
num_subjects = 6;
num_states = zeros(num_subjects, size(trialset,2));
conflict_all = []; 
total_states_all = []; 
subject_ids = [];

% Create a color gradient from black to light grey for each subject.
colors = zeros(num_subjects, 3);
for n = 1:num_subjects
    c = (n-1) * (0.8/(num_subjects-1));  % subject 1 is [0,0,0] and subject 6 is [0.8,0.8,0.8]
    colors(n,:) = [c, c, c];
end

% Set the jitter magnitude (adjust as needed)
jitterMagnitude = 0.01;

for n = 1:num_subjects
    clear trialset;
    unique_trial_types = unique(subject(n).trial_type_trial);
    for i = 1:length(unique_trial_types)
        trialset{i} = find(subject(n).trial_type_trial == unique_trial_types(i) & subject(n).decision < 2);
    end

    state_metadata_name = strcat(['state_metadata_' num2str(n) '.mat']);
    load(state_metadata_name)

    for a = 1:length(trialset)
        filter = trialset{a};
        total_states = 0;
        for i = 1:length(filter)
            entry_times_trial = entry_times(trial_identity == filter(i));
            total_states = total_states + length(entry_times_trial);
        end
        % Average the number of states across trials for this trial type
        num_states(n,a) = total_states / length(filter);
    end

    % Use the original conflict values and add jitter for visualization.
    x_original = subject(n).conflict_trial_type;
    x_jittered = x_original + jitterMagnitude * randn(size(x_original));
    
    % Scatter plot with 40% opacity
    scatter(x_jittered, num_states(n, unique_trial_types), 36, colors(n,:), ...
        'filled', 'MarkerFaceAlpha', 0.4)
    hold on
    [h_ind(n), p_ind(n)] = corr(x_original', num_states(n, unique_trial_types)');
    
    subject_ids = [subject_ids, ones(1, length(unique_trial_types)) * n];
    conflict_all = [conflict_all, x_original];
    total_states_all = [total_states_all, num_states(n, unique_trial_types)];
    
    % Compute the best-fit line using the original (non-jittered) conflict values.
    p_coeff = polyfit(x_original, num_states(n, unique_trial_types), 1);
    x_fit = linspace(min(x_original), max(x_original), 100);
    y_fit = polyval(p_coeff, x_fit);
    
    % Plot the best-fit line using the same subject-specific color.
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', colors(n,:));
end

xlabel('Conflict')
ylabel('Number of States Visited')
title('Best Fit Lines for Each Subject')
%% stats- linear mixed effects model - does conflict modulate number of states visited?

tbl = table(subject_ids', conflict_all', total_states_all',...
    'VariableNames', {'Subject', 'Conflict', 'TotalStates'});

tbl.Subject = categorical(tbl.Subject);

% Define a mixed-effects model:
% The fixed effect is 'Conflict', and we include a random intercept for each subject.
lme = fitlme(tbl, 'TotalStates ~ Conflict + Subject');

% View the model results:
disp(lme)
