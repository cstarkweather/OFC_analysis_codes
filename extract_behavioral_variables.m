% Define the directories
folders = {'EC304', 'EC288', 'PR05', 'PR06','BJH058','DP01'};

for n=1:length(folders)
processFile(folders{n});
end

%%
% Function to process each file
function processFile(filePath)
    cd(filePath)
    behavior_file = dir('*.txt');
    behavior_file.name
    % Import the data without treating the first row as a header
    opts = detectImportOptions(behavior_file.name, 'NumHeaderLines', 0);
    tableData = readtable(behavior_file.name, opts);


    mat = NaN(5, 8);
    reward_sizes = 0:7;
    punishment_sizes = 0:4;


    % Extract VAS from the actual first row
    VAS = tableData{1, :};

    % Skip the first row to process the rest of the data
    tableData(1, :) = [];

    % Use since not using triggers
    trial_type = tableData.Var2;
    decision = tableData.Var3 > 0;
    reaction_time = tableData.Var7 - tableData.Var5;

    % Figure out whether 192 vs 220 trial version of task (more granular trials in 192 trial version)
    if length(unique(trial_type)) < 12
        trialID = [1 2 3 4 5 6 7 8 9 10 11];
        p = [0 2 4 0 2 4 0 2 4 2 4];
        r = [2 2 2 4 4 4 7 7 7 0 0];
    else
        trialID = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
        p = [0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 2 4];
        r = [2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 7 7 7 7 7 0 0];
    end

    clear teleported
    clear teleported_all
    expanded_predictors = zeros(length(decision), 2);

    for i = 1:length(decision)
        expanded_predictors(i, :) = [p(trial_type(i)) r(trial_type(i))];
    end
     pun = expanded_predictors(:,1);rew = expanded_predictors(:,2);
    for i = 1:length(trialID)
        indices = find(trial_type == i);
        approached(i) = sum(decision(indices) == 1) / length(indices);
        approached_all{i} = decision(indices) == 1;
    end

    % Fill the matrix with your data
    for i = 1:5 % cycle through different punishment sizes
        for j = 1:8
            index = find(p == punishment_sizes(i) & r == reward_sizes(j));
            if ~isempty(index)
                mat(i, j) = approached(index);
                mat_all{i, j} = approached_all{i};
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
[b, dev, stats] = glmfit([X_flat Y_flat], P_flat, 'binomial', 'link', 'logit', 'constant', 'off');

% Display the coefficients

%a1 = bLasso(1);  % Coefficient for reward_sizes (x)
%a2 = bLasso(2);  % Coefficient for punishment_sizes (y)

    reward_coeff = b(1);
    punish_coeff = b(2);
    intercept_coeff = 0;%interceptLasso; % this is also val_av

    % Compute entropy for each item in P_flat
for i = 1:length(P_flat)
    prob = P_flat(i);
    
    % Handle the special cases where p is exactly 0 or 1 to avoid log(0)
    if prob == 0 || prob == 1
        conflict((find(r == X_flat(i) & p == Y_flat(i)))) = 0;  % Entropy is 0 when the outcome is certain
    else
        % Compute the entropy using the formula
        conflict((find(r == X_flat(i) & p == Y_flat(i)))) = -prob * log2(prob) - (1 - prob) * log2(1 - prob);
    end
end

p_approach_trial = zeros(length(trialID),1);
for i=1:length(trial_type)
    p_approach_trial(i)= approached(trialID==trial_type(i));
end

conflict_trial = conflict;
reward_trial = rew;
punish_trial = pun;
value_trial = reward_coeff * reward_trial + punish_coeff * punish_trial;
trial_type_trial = trial_type;

reward_trial_type = r;
punishment_trial_type = p;
conflict_trial_type = conflict;
trial_type_trial_type = trialID;
value_trial_type = reward_coeff * reward_trial_type + punish_coeff * punishment_trial_type;
p_approach_trial_type = approached;

save('behavior.mat','VAS','reaction_time','reward_coeff','punish_coeff','p_approach_trial','conflict_trial','reward_trial','punish_trial','value_trial','trial_type_trial',...
    'reward_trial_type','punishment_trial_type','conflict_trial_type','trial_type_trial_type','value_trial_type','p_approach_trial_type')

% Return to the parent directory
cd('..');

end