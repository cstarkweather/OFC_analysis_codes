
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
                if n ==3 && array ==1 %plotting example matrix
                    figure(1)
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

            else
            end


        end
    end

fprintf('Total number of negative correlated pairs:')
size(neg_corrs,1)

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

%compare distribution across clusters for the more laterally located
%electrode in each negatively correlated pair
observed = lateral_in_clust * size(neg_corrs,1);      % Arbitrary total N=100
expected = lateral_in_clust_null * size(lateral,1);   % Same arbitrary total
fprintf('comparison of distribution across clusters for more laterally located electrode in each negatively correlated pair')
[h, p, stats] = chi2gof([1,2,3], 'Freq', observed, 'Expected', expected, 'Emin', 1)
fprintf('Chi-square(%d) = %.2f, p = %.2e\n', stats.df, stats.chi2stat, p);

%compare distribution across clusters for the more medially located
%electrode in each negatively correlated pair

observed = medial_in_clust * size(neg_corrs,1);      % Arbitrary total N=100
expected = medial_in_clust_null * size(medial,1);    % Same arbitrary total
fprintf('comparison of distribution across clusters for more laterally located electrode in each negatively correlated pair')
[h, p, stats] = chi2gof([1,2,3], 'Freq', observed, 'Expected', expected, 'Emin', 1)
fprintf('Chi-square(%d) = %.2f, p = %.2e\n', stats.df, stats.chi2stat, p);

%% heatmap of spatial distribution of the more medially located electrode in each negatively correlated pair
figure(2)
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
X_coords = neg_corrs(:, 4);
Y_coords = neg_corrs(:, 6);

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
%% heatmap of spatial distribution of the more laterally located electrode in each negatively correlated pair
figure(3)
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

%% compare null distribution of the medial/lateral distance between the medial and lateral electrodes with the actual distribution
null_lat = null_corrs(:,5) - null_corrs(:,4);
neg_lat = neg_corrs(:,5) - neg_corrs(:,4);

[ a b c d] = anovan([null_lat;neg_lat],[zeros(length(null_lat),1);ones(length(neg_lat),1)]);
 
%% compare null distribution of the anterior/posterior distance between the medial and lateral electrodes with the actual distribution
null_ap = null_corrs(:,7) - null_corrs(:,6);
neg_ap = neg_corrs(:,7) - neg_corrs(:,6);

anovan([null_ap;neg_ap],[zeros(length(null_ap),1);ones(length(neg_ap),1)]);


