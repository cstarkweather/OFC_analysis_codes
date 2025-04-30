
%% First organize electrodes that are all on on hemisphere of each subject's brain
% SUBJECT 1
subject(1).array=zeros(length(subject(1).electrode),1);
for i=1:5
subject(1).array(i) = 1;
end
for i=6:10
subject(1).array(i) = 2;
end

% SUBJECT 2
subject(2).array=zeros(length(subject(2).electrode),1);
for i=3:9
subject(2).array(i) = 1;
end
for i=11:17
subject(2).array(i) = 2;
end

% SUBJECT 3
subject(3).array=zeros(length(subject(3).electrode),1);
for i=4:15
subject(3).array(i) = 1;
end
for i=16:24
subject(3).array(i) = 1;
end
for i=28:39
subject(3).array(i) = 2;
end

% SUBJECT 4
subject(4).array=zeros(length(subject(4).electrode),1);
for i=3:14
subject(4).array(i) = 1;
end
for i=16:20
subject(4).array(i) = 2;
end

% SUBJECT 5
subject(5).array=zeros(length(subject(5).electrode),1);
for i=8:17
subject(5).array(i) = 1;
end

% SUBJECT 6
subject(6).array=zeros(length(subject(6).electrode),1);
for i=3:10
subject(6).array(i) = 1;
end
for i=13:20
subject(6).array(i) = 2;
end

%% 

bound1 = 4700; bound2 = 5000;
total_comparisons = 0;
neg_corrs=[];noise_correlation_store = [];null_corrs=[];


    neg_corrs=[];noise_correlation_store = [];
    for array = 1:2
        % Initialize noise_electrode cell array
        for n = 1:6
            clear I;
            clear coordinate_med;clear coordinate_trans;clear noise_electrode;
            electrodes = find(subject(n).array==array);
            clear noise_electrode;

            which_trials1 = find(subject(n).decision == 0);
            which_trials2 = find(subject(n).decision == 1);

            % Compute noise correlation matrix and p-values
            num_electrodes = length(electrodes);

            counter = 1;

            for i = 1:num_electrodes
                e = electrodes(i);

            % with a baseline
            baseline1 = sum(subject(n).electrode(e).trigger(1).high_gamma_mat(which_trials1,4801:5000),2)/200;
            baseline2 = sum(subject(n).electrode(e).trigger(1).high_gamma_mat(which_trials2,4801:5000),2)/200;

            avg_1 = sum((subject(n).electrode(e).trigger(2).high_gamma_mat(which_trials1,bound1:bound2)-baseline1)./baseline1) / length(which_trials1);
            avg_2 = sum((subject(n).electrode(e).trigger(2).high_gamma_mat(which_trials2,bound1:bound2)-baseline2)./baseline2) / length(which_trials2);

            noise1 = (subject(n).electrode(e).trigger(2).high_gamma_mat(which_trials1,bound1:bound2)-baseline1)./(baseline1) - avg_1;
            noise2 = (subject(n).electrode(e).trigger(2).high_gamma_mat(which_trials2,bound1:bound2)-baseline2)./(baseline2) - avg_2;

            noise_electrode{counter} = [noise1;noise2]; % Combine the noise matrices
            coordinate_med(counter) = subject(n).electrode(e).med_coor;
            coordinate_trans(counter) = subject(n).electrode(e).trans_coor;
            counter = counter+1;

            end

            noise_correlation_matrix = zeros(counter-1, counter-1);
            pval_matrix = zeros(counter-1, counter-1);

            if exist('coordinate_med')
                [B,I] = sort(coordinate_med);

                for i = 1:counter-1
                    for j = i:counter-1
                        if i ==j | abs(i - j)<1
                            noise_correlation_matrix(I(i), I(j)) = 1;
                            pval_matrix(I(i), I(j)) = 1;
                        else

                            % Compute Pearson correlation and p-value
                            [corr_value, p_value] = calculate_corr(noise_electrode{I(i)}, noise_electrode{I(j)});

                            % Store in matrices
                            noise_correlation_matrix(i,j) = corr_value;
                            noise_correlation_matrix(j,i) = corr_value; % Symmetric matrix
                            if corr_value < 0
                                neg_corrs = [neg_corrs; n electrodes(I(i)) electrodes(I(j)) B(i) B(j) coordinate_trans(I(i)) coordinate_trans(I(j)) p_value corr_value];
                            else

                            end

                            null_corrs = [null_corrs; n electrodes(I(i)) electrodes(I(j)) B(i) B(j) coordinate_trans(I(i)) coordinate_trans(I(j)) corr_value];
                            total_comparisons = total_comparisons + 1;

                            pval_matrix(i,j) = p_value;
                            pval_matrix(j,i) = p_value; % Symmetric matrix
                        end
                    end
                end

                % Convert the matrix into a single column vector
                vector_column1 = noise_correlation_matrix(:);

                % Create a second column with the value of n repeated
                vector_column2 = repmat(n, length(vector_column1), 1);

                % Combine into a two-column matrix
                result_matrix = [vector_column1, vector_column2];

                noise_correlation_store = [noise_correlation_store;result_matrix];

                noise_correlation_matrices{n} = noise_correlation_matrix;
                pval_matrices{n} = log10(pval_matrix);

                % Subplot 1: Noise correlation matrix
                subplot(3,6, n+(array-1)*6);
                imagesc(noise_correlation_matrices{n}); % Display correlation matrix
                colorbar;
                colormap gray
                caxis([-0.15 0])
                all_electrodes = electrodes(I);
                num_electrodes = length(all_electrodes);
                set(gca, 'XTick', 1:num_electrodes, 'XTickLabel', all_electrodes);
                set(gca, 'YTick', 1:num_electrodes, 'YTickLabel', all_electrodes);
                axis square;
                box off
                set(gca, 'FontSize', 8); % Reduce font size

                hold on;

            else
            end


        end
    end


size(neg_corrs)
%% plot medial/lateral coordinates of each pair
subplot(2,1,1)
histogram(neg_corrs(1:end,4),'binwidth',2)
xlim([-20 50])
subplot(2,1,2)
histogram(neg_corrs(1:end,5),'binwidth',2)
xlim([-20 50])
%% plot anterior/posterior coordinates of each pair
subplot(2,1,1)
histogram(neg_corrs(:,6),10)
xlim([-15 50])
subplot(2,1,2)
histogram(neg_corrs(:,7),10)
xlim([-15 50])

%% number of highly significant (p < 0.001) negative correlations
size(neg_corrs(neg_corrs(:,8) < 0.001,:))
unique(neg_corrs(neg_corrs(:,8) < 0.001,1))

%% what percent of Cluster 1 is one of these more medial electrodes?

for i = 1
    % Get electrode and subject IDs for the current cluster
    cluster_electrodes = flagged_electrode_ID(clustlabel == i);
    cluster_subjects = flagged_subject_ID(clustlabel == i);

    % Combine electrodes and subjects into pairs for the cluster
    cluster_pairs = [cluster_subjects, cluster_electrodes'];

    % Combine electrodes and subjects into pairs for negative coding
    medial = [neg_corrs(:,1), neg_corrs(:,2)];

    % Find the indices of negative coding electrodes that are in the cluster
    medial_in_clust = ismember(cluster_pairs,medial, 'rows');

    % Number of negative coding electrodes in cluster
    medial_located_in_clust = sum(medial_in_clust);

    % Proportion of electrodes in cluster that correspond to the medially
    % located, negatively correlated group
    proportion_medial_in_clust = medial_located_in_clust/ sum(clustlabel == i)

    % Proportion of electrodes in cluster that correspond to the medially
    % located, negatively correlated group
    proportion_medial_in_clust = medial_located_in_clust/ sum(clustlabel == i)

end

%% what percent of the more medial paired lies in Cluster 1?

for i = 1:3
    % Get electrode and subject IDs for the current cluster
    cluster_electrodes = flagged_electrode_ID(clustlabel == i);
    cluster_subjects = flagged_subject_ID(clustlabel == i);

    % Combine electrodes and subjects into pairs for the cluster
    cluster_pairs = [cluster_subjects, cluster_electrodes'];

    % Combine electrodes and subjects into pairs for negative coding
    medial = [neg_corrs(:,1), neg_corrs(:,2)];
    sum(ismember(medial,cluster_pairs, 'rows'))
    % Find the indices of anticorrelated, medially located electrodes that are in the cluster
    medial_in_clust(i) = sum(ismember(medial,cluster_pairs, 'rows'))/size(medial,1);
    medial_in_clust_null(i) = sum(clustlabel==i)/length(clustlabel);

    % Combine electrodes and subjects into pairs for negative coding
    lateral = [neg_corrs(:,1), neg_corrs(:,3)];

    % Find the indices of anticorrelated, medially located electrodes that are in the cluster
    lateral_in_clust(i) = sum(ismember(lateral,cluster_pairs, 'rows'))/size(lateral,1);
    lateral_in_clust_null(i) = sum(clustlabel==i)/length(clustlabel);

end

observed = lateral_in_clust * size(neg_corrs,1)      % Arbitrary total N=100
expected = lateral_in_clust_null * size(lateral,1) % Same arbitrary total

[h, p, stats] = chi2gof([1,2,3], 'Freq', observed, 'Expected', expected, 'Emin', 1)
fprintf('Chi-square(%d) = %.2f, p = %.2e\n', stats.df, stats.chi2stat, p);


observed = medial_in_clust * size(neg_corrs,1)      % Arbitrary total N=100
expected = medial_in_clust_null * size(medial,1) % Same arbitrary total

[h, p, stats] = chi2gof([1,2,3], 'Freq', observed, 'Expected', expected, 'Emin', 1)
fprintf('Chi-square(%d) = %.2f, p = %.2e\n', stats.df, stats.chi2stat, p);




%% heatmap
% Parameters
window_size = 10;  % Size of the sliding window in mm
overlap = 5;       % Overlap between windows in mm

% Extract X and Y coordinates
X_coords = neg_corrs(:, 4);
Y_coords = neg_corrs(:, 6);

% Determine the min and max coordinates for the grid
x_min = floor(min(X_coords)) - window_size;
x_max = ceil(max(X_coords)) + window_size;
y_min = floor(min(Y_coords)) - window_size;
y_max = ceil(max(Y_coords)) + window_size;

% Generate sliding window centers with overlap
x_centers = x_min:overlap:x_max;
y_centers = y_min:overlap:y_max;

% Initialize matrix to store electrode counts in each bin
electrode_counts = zeros(length(y_centers), length(x_centers));

% Extract X and Y coordinates
X_coords = neg_corrs(:, 5);
Y_coords = neg_corrs(:, 7);

% Loop over sliding window centers to calculate the electrode count in each bin
for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        % Define the window range
        x_start = x_centers(x_idx) - window_size/2;
        x_end = x_centers(x_idx) + window_size/2;
        y_start = y_centers(y_idx) - window_size/2;
        y_end = y_centers(y_idx) + window_size/2;

        % Count electrodes within the current window
        in_window = X_coords >= x_start & X_coords < x_end & ...
                    Y_coords >= y_start & Y_coords < y_end;
        electrode_counts(y_idx, x_idx) = sum(in_window);
    end
end

% Set bins with no electrodes to NaN for white background
plot_data = electrode_counts;
plot_data(electrode_counts == 0) = NaN;  % Set empty bins to NaN

% Plot the heatmap
figure;
h = imagesc(x_centers, y_centers, plot_data);  % Plot the data
colormap([0 0 0 ]);  % Set hot pink as the colormap
%colorbar off;  % Turn off the colorbar since color doesn't vary

% Set transparency based on the number of electrodes in each bin
alpha_data = electrode_counts ./ max(max(electrode_counts)); % Normalize alpha
set(h, 'AlphaData', alpha_data);

% Set axis limits
xlim([-25, 45]);
ylim([-25, 30]);

% Ensure y-axis is not flipped
set(gca, 'YDir', 'normal');

% Add labels and title
title('Electrode Density Heatmap', 'FontSize', 14);
xlabel('X Coordinate (mm)', 'Fontsize', 14);
ylabel('Y Coordinate (mm)', 'Fontsize', 14);
set(gca, 'TickDir', 'out');
box off;
hold on;

%% euclidean distances
figure(2)
% Parameters
window_size = 10;  % Size of the sliding window in mm
overlap = 2.5;       % Overlap between windows in mm

% Extract X and Y coordinates for both electrodes in the pair
X1 = neg_corrs(:, 4);  % X coordinate of the first electrode
Y1 = neg_corrs(:, 6);  % Y coordinate of the first electrode
X2 = neg_corrs(:, 5);  % X coordinate of the second electrode
Y2 = neg_corrs(:, 7);  % Y coordinate of the second electrode

% Compute signed distances for each pair
delta_X = X2 - X1;  % Signed difference in X coordinates
delta_Y = Y2 - Y1;  % Signed difference in Y coordinates

% Determine the min and max coordinates for the grid
x_min = floor(min(delta_X)) - window_size;
x_max = ceil(max(delta_X)) + window_size;
y_min = floor(min(delta_Y)) - window_size;
y_max = ceil(max(delta_Y)) + window_size;

% Generate sliding window centers with overlap
x_centers = x_min:overlap:x_max;
y_centers = y_min:overlap:y_max;

% Initialize matrix to store counts of pairs in each bin
signed_distance_counts = zeros(length(y_centers), length(x_centers));

% Loop over sliding window centers to calculate counts of signed distances in each bin
for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        % Define the window range
        x_start = x_centers(x_idx) - window_size/2;
        x_end = x_centers(x_idx) + window_size/2;
        y_start = y_centers(y_idx) - window_size/2;
        y_end = y_centers(y_idx) + window_size/2;

        % Count pairs within the current window
        in_window = delta_X >= x_start & delta_X < x_end & ...
                    delta_Y >= y_start & delta_Y < y_end;
        signed_distance_counts(y_idx, x_idx) = sum(in_window);
    end
end

% Plot the heatmap
figure;
h = imagesc(x_centers, y_centers, signed_distance_counts);  % Plot the data
colormap('cool');  % Use a hot colormap for visualization
colorbar;
xlabel('ΔX (mm)', 'FontSize', 12);
ylabel('ΔY (mm)', 'FontSize', 12);
title('Heatmap of Signed Distances Between Electrode Pairs', 'FontSize', 14);

% Ensure y-axis is not flipped
set(gca, 'YDir', 'normal');

% Add transparency to reflect missing data
alpha_data = signed_distance_counts ./ max(max(signed_distance_counts)); % Normalize alpha
set(h, 'AlphaData', alpha_data);

% Set axis properties
set(gca, 'TickDir', 'out');
box off;
%%
colormap(flipud(gray(1000))); % 24 levels from black to white (0–23)
colorbar;
caxis([23 0]);      % Set color axis limits
%% euclidean distances - null corrs
% Parameters
window_size = 10;  % Size of the sliding window in mm
overlap = 2.5;       % Overlap between windows in mm

% Extract X and Y coordinates for both electrodes in the pair
X1 = null_corrs(:, 4);  % X coordinate of the first electrode
Y1 = null_corrs(:, 6);  % Y coordinate of the first electrode
X2 = null_corrs(:, 5);  % X coordinate of the second electrode
Y2 = null_corrs(:, 7);  % Y coordinate of the second electrode

% Compute signed distances for each pair
delta_X = X2 - X1;  % Signed difference in X coordinates
delta_Y = Y2 - Y1;  % Signed difference in Y coordinates

% Determine the min and max coordinates for the grid
x_min = floor(min(delta_X)) - window_size;
x_max = ceil(max(delta_X)) + window_size;
y_min = floor(min(delta_Y)) - window_size;
y_max = ceil(max(delta_Y)) + window_size;

% Generate sliding window centers with overlap
x_centers = x_min:overlap:x_max;
y_centers = y_min:overlap:y_max;

% Initialize matrix to store counts of pairs in each bin
signed_distance_counts = zeros(length(y_centers), length(x_centers));

% Loop over sliding window centers to calculate counts of signed distances in each bin
for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        % Define the window range
        x_start = x_centers(x_idx) - window_size/2;
        x_end = x_centers(x_idx) + window_size/2;
        y_start = y_centers(y_idx) - window_size/2;
        y_end = y_centers(y_idx) + window_size/2;

        % Count pairs within the current window
        in_window = delta_X >= x_start & delta_X < x_end & ...
                    delta_Y >= y_start & delta_Y < y_end;
        signed_distance_counts(y_idx, x_idx) = sum(in_window);
    end
end

% Plot the heatmap
figure;
h = imagesc(x_centers, y_centers, signed_distance_counts);  % Plot the data
colormap('cool');  % Use a hot colormap for visualization
colorbar;
xlabel('ΔX (mm)', 'FontSize', 12);
ylabel('ΔY (mm)', 'FontSize', 12);
title('Heatmap of Signed Distances Between Electrode Pairs', 'FontSize', 14);

% Ensure y-axis is not flipped
set(gca, 'YDir', 'normal');

% Add transparency to reflect missing data
alpha_data_null = signed_distance_counts ./ max(max(signed_distance_counts)); % Normalize alpha
set(h, 'AlphaData', alpha_data_null);

% Set axis properties
set(gca, 'TickDir', 'out');
box off;
%%
binwidth = 5;
subplot(4,1,1)
histogram(null_corrs(:,7) - null_corrs(:,6),'Binwidth',binwidth)
xlim([-50 50])
subplot(4,1,2)
histogram(neg_corrs(:,7) - neg_corrs(:,6),'Binwidth',binwidth)
xlim([-50 50])
subplot(4,1,3)
histogram(null_corrs(:,5) - null_corrs(:,4),'Binwidth',binwidth)
xlim([-50 50])
subplot(4,1,4)
histogram(neg_corrs(:,5) - neg_corrs(:,4),'Binwidth',binwidth)
xlim([-50 50])
%%
null_lat = null_corrs(:,5) - null_corrs(:,4);
neg_lat = neg_corrs(:,5) - neg_corrs(:,4);

[ a b c d] = anovan([null_lat;neg_lat],[zeros(length(null_lat),1);ones(length(neg_lat),1)])
sum(neg_lat)/length(neg_lat)
sum(null_lat)/length(null_lat) 
%%
null_ap = null_corrs(:,7) - null_corrs(:,6);
neg_ap = neg_corrs(:,7) - neg_corrs(:,6);

anovan([null_ap;neg_ap],[zeros(length(null_ap),1);ones(length(neg_ap),1)])
sum(neg_ap)/length(neg_ap)
sum(null_ap)/length(null_ap) 
%%
neg_corrs_1

%% subplots of delta X and delta Y
subplot(2,1,1)
% Plot individual points for delta_X
for i = 1:size(X1, 1)  % Use size(X1, 1) to loop over rows
    plot([rand*0.1 + 1], [delta_X(i)], 'k.', 'Markersize', 10)
    hold on
end

% Compute and plot mean and standard deviation for delta_X
mean_delta_X = mean(delta_X);
std_delta_X = std(delta_X);
errorbar(1, mean_delta_X, std_delta_X, 'r', 'LineWidth', 1.5, 'CapSize', 0); % Mean and std dev
xlim([0.5 1.5])
set(gca, 'TickDir', 'out');
box off;
title('Delta X with Mean and Standard Deviation');
ylabel('Delta X');

subplot(2,1,2)
% Plot individual points for delta_Y
for i = 1:size(X1, 1)  % Use size(X1, 1) to loop over rows
    plot([rand*0.1 + 1], [delta_Y(i)], 'k.', 'Markersize', 10)
    hold on
end

% Compute and plot mean and standard deviation for delta_Y
mean_delta_Y = mean(delta_Y);
std_delta_Y = std(delta_Y);
errorbar(1, mean_delta_Y, std_delta_Y, 'r', 'LineWidth', 1.5, 'CapSize', 0); % Mean and std dev
xlim([0.5 1.5])
set(gca, 'TickDir', 'out');
box off;
title('Delta Y with Mean and Standard Deviation');
ylabel('Delta Y');

%%
subplot(2,1,2)
plot([rand(size(Y1))*0.1 + 1;rand(size(Y2))*0.1 + 2],[Y1;Y2],'k.','Markersize',10)
xlim([0.5 2.5])
set(gca, 'TickDir', 'out');
box off;
%% colorbar for Subject 3
% Extract the range for the colorbar
min_value = subject(3).electrode(4).med_coor;
max_value = subject(3).electrode(15).med_coor;

heatmap_data = zeros(12,1);
for i=1:12
    heatmap_data(i) = subject(3).electrode(i+3).med_coor;
end

% Plot the heatmap
figure;
imagesc(heatmap_data);  % Plot the data
colormap('spring');        % Apply the jet colormap
colorbar;

% Set the color axis limits to span the desired range
caxis([min_value, max_value]);
%% colorbar for Subject 5
% Extract the range for the colorbar
min_value = subject(3).electrode(8).med_coor;
max_value = subject(3).electrode(17).med_coor;

heatmap_data = zeros(10,1);
for i=1:10
    heatmap_data(i) = subject(5).electrode(i+7).med_coor;
end

% Plot the heatmap
figure;
imagesc(heatmap_data);  % Plot the data
colormap('spring');        % Apply the jet colormap
colorbar;

% Set the color axis limits to span the desired range
caxis([min_value, max_value]);


%% look at group PSTH

t = 3;num_conditions = 2;

signal_bounds = [3000 7000];
baseline_bounds = [4800 5000];
window = signal_bounds(2) - signal_bounds(1) + 1;
smoothing_kernel = 1;
x2 = round(signal_bounds(1):smoothing_kernel:signal_bounds(2));
colorflag = 0;

for j = 1:num_conditions

    for_psth = zeros(length(clust_indices), window);

    for i = 1:size(neg_corrs,1)

        n = neg_corrs(i,1);
        e = neg_corrs(i,2);

        trialset = [{subject(n).decision==0} {subject(n).decision==1}];

        trialIDs=trialset{j};

        % 'temp' gathers what the electrode does for all trials
        temp = [subject(n).electrode(e).trigger(t).high_gamma_mat(trialIDs, signal_bounds(1):signal_bounds(2))];

        % baseline
        temp_baseline = [subject(n).electrode(e).trigger(1).high_gamma_mat(trialIDs, baseline_bounds(1):baseline_bounds(2))];
        baseline = (sum(sum(temp_baseline))/(size(temp_baseline,1)* size(temp_baseline,2)));

        for_psth(i, :) = (sum(temp,1)/size(temp,1)-baseline)/baseline;

    end

    % Remove rows that contain any NaN values
    for_psth = for_psth(~any(isnan(for_psth), 2), :);


    y = smooth(sum(for_psth) / size(for_psth, 1), smoothing_kernel);
    dy = zeros(length(y),1);
    for c = 1:length(x2) - 1
        dy((c - 1) * smoothing_kernel + 1:c * smoothing_kernel) = std(sum(for_psth(:, (c - 1) * smoothing_kernel + 1:c * smoothing_kernel),2)/smoothing_kernel)/sqrt(size(for_psth,1));
    end

    for_anova{j} = for_psth;

    smooth_y = y(1:smoothing_kernel:end+1-smoothing_kernel)';
    smoothed_dy = smooth(dy,smoothing_kernel);
    if colorflag ==0
        fill([(signal_bounds(1):signal_bounds(2))-5000, fliplr((signal_bounds(1):signal_bounds(2))-5000)], [(y - smoothed_dy)'  fliplr((y + smoothed_dy)')],[1 - (j - 1) / (length(trialset) - 0.9), (j - 1) / (length(trialset) - 0.9), 0.15], 'linestyle', 'none','FaceAlpha', 0.6);
        hold on
        plot((signal_bounds(1):signal_bounds(2))'-5000, y, 'Linewidth',2,'Color',[0.7*(1 - (j - 1) / (length(trialset) - 1)) 0.5*(j - 1) / (length(trialset) - 1)  0 1])
        hold on
    else

    fill([(signal_bounds(1):signal_bounds(2))-5000, fliplr((signal_bounds(1):signal_bounds(2))-5000)], [(y - smoothed_dy)'  fliplr((y + smoothed_dy)')], [1 - (j-1) / (length(trialset)-1), (j -1)/ (length(trialset)-1), 1], 'linestyle', 'none','FaceAlpha', 0.6);
    
    hold on
    plot((signal_bounds(1):signal_bounds(2))'-5000, y, 'k','Linewidth',1)
    hold on
    end

end


% Adjust plot appearance
set(gca, 'TickDir', 'out');
% Get the current x-axis limits
x_limits = xlim;

% Define tick marks, ensuring that 0 is included
x_ticks = linspace(x_limits(1), x_limits(2), 2);  % Create 5 evenly spaced ticks
x_ticks = unique([x_ticks, 0]);  % Ensure 0 is included

% Set the x-axis ticks
set(gca, 'XTick', x_ticks);


% Get the current y-axis limits
y_limits = ylim;

% Define tick marks, ensuring that 0 is included
y_ticks = linspace(y_limits(1), y_limits(2), 2);  % Create 5 evenly spaced ticks
y_ticks = unique([y_ticks, 0]);  % Ensure 0 is included
y_ticks = unique([y_ticks, 0.1]);  % Ensure 0 is included

% Set the y-axis ticks
set(gca, 'YTick', y_ticks);


box off;
set(gcf, 'Color', 'w');
xlabel('time from decision (milliseconds)', 'FontSize', 20);
ylabel('\Delta High Gamma Amplitude', 'Interpreter', 'tex', 'FontSize', 20);

xline(0,'k--','Linewidth',2)