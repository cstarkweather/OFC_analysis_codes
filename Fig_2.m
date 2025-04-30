
%% build roc_tc
folders = [{'EC304'} {'EC288'} {'PR05'} {'PR06'} {'BJH058'} {'DP01'}];

% Initialize variables for concatenation
all_roc_tc = [];
olf_coordinates = [];
med_coordinates = [];
a_p_coordinates=[];
all_subjects = {};
electrode_ID =[];
subject_ID=[];
side=[];brain_region=[];
lat_coordinates = [];

% Parameters for ROC calculation
bin_width = 50; % User-specified bin width in milliseconds
t = 2;

% Loop through each subject
for n = 1:length(subject)
    % Loop through each electrode
    for e = 1:length(subject(n).electrode)

        signal = subject(n).electrode(e).trigger(t).high_gamma_mat(:,4000:6000);

        trialIDs1 = find(subject(n).decision==1); % User-specified trial IDs, 1 for approach trials
        binned_data1 = bin_data_columns(bin_width, trialIDs1, signal);

        signal_baseline_1 = subject(n).electrode(e).trigger(1).high_gamma_mat(trialIDs1, 4800:5000); % Unwrap along rows

        basedistr = bin_data_columns(bin_width, 1:size(signal_baseline_1,1), signal_baseline_1);
        roc_tc1 = zeros(1,size(binned_data1, 2));
        % Compute ROC values for each row of binned_data and basedistr
        for j = 1:size(binned_data1, 2)
            distr = binned_data1(:, j);
            roc_tc1(j) = auROC(distr, basedistr, 1000);
        end


        all_roc_tc=[all_roc_tc;roc_tc1];
        side = [side;strcmp(subject(n).electrode(e).brain_region(1),'L')];
        brain_region = [brain_region;strcmp(subject(n).electrode(e).brain_region,'RSGC')];
        olf_coordinates = [olf_coordinates; subject(n).electrode(e).olf_coor];  
        med_coordinates = [med_coordinates; subject(n).electrode(e).med_coor]; 
        lat_coordinates = [lat_coordinates; subject(n).electrode(e).lat_coor]; 
        a_p_coordinates = [a_p_coordinates; subject(n).electrode(e).trans_coor];  
        all_subjects = [all_subjects; folders{n}]; 
        subject_ID = [subject_ID; n]; 
        electrode_ID  = [electrode_ID e];
    end
end

%% generate flags for excluding any electrodes

% Exclude certain indices

exclude_indices = [olf_coordinates < -1];

%for excluding arrays with majority posterior electrodes
%exclude_indices = [olf_coordinates < -1 | (subject_ID==6) | (subject_ID==4 & side==1) | (subject_ID==3 & brain_region==1)];

% Create new variables to store the flagged data
flagged_roc_tc = all_roc_tc;
flagged_coordinates = med_coordinates;
flagged_trans_coordinates = a_p_coordinates;
flagged_lat_coordinates = lat_coordinates;
flagged_subjects = all_subjects;
flagged_subject_ID = subject_ID;
flagged_electrode_ID = electrode_ID;
flagged_side = side;

% Remove the excluded indices from flagged_roc_tc and flagged_coordinates
flagged_roc_tc(exclude_indices, :) = [];
flagged_coordinates(exclude_indices) = [];
flagged_lat_coordinates(exclude_indices) = [];
flagged_trans_coordinates(exclude_indices) = [];

flagged_subjects(exclude_indices) = [];
flagged_subject_ID(exclude_indices) = [];
flagged_electrode_ID(exclude_indices) = [];
flagged_side(exclude_indices) = [];




%% Determine plot order - if using K means clustering
% PCA implementation: pca().
% [COEFF, SCORE] = PCA(X) returns the principal component scores,
%     i.e., the projections of X in the principal component space.  Rows
%     of SCORE correspond to observations, columns to components.
% Coeff is a column matrix containing ordered principal components:
%   e.g., coeff(:,1) is the first PC.
% Scores are the projections of each datapoint into PC-space.  e.g.,
%   score(:,1) is the first projection
% Latency are  eigenvalues
sorting_attribute = flagged_roc_tc;pc_roc=NaN;
[pc_roc, proj_roc, eigv_roc,~,explained] = ...
    pca([sorting_attribute(:,1:21)]);

plot(pc_roc(:,1));
legend({'PC1'});
eigv_roc(1:3)  %display the values of the first three eigenvalues

% Compute Kmeans cluster, relabel and plot the clustering result.
nclust = 3;
data = proj_roc(:,1); %projections of data onto PCs
seeds = kmeansinit(data,nclust);
clustlabel = kmeans(data,nclust, 'Start',seeds);


box off;
set(gca, 'TickDir', 'out');  % Set tick marks to face out

[~, plotorder] = sort(clustlabel);

sorted_roc_tc = sorting_attribute(plotorder, :);
sorted_subjects = flagged_subjects(plotorder);
sorted_coordinates = flagged_coordinates(plotorder);
sorted_subject_ID = flagged_subject_ID(plotorder);
sorted_electrode_ID = flagged_electrode_ID(plotorder);

% %% if want to plot pc scores
%imagesc(proj_roc(plotorder,1));colormap bone;box off;
%set(gca, 'TickDir', 'out');  % Set tick marks to face out
%colorbar
% %% % variance explained
% plot(explained)
% box off
% set(gca, 'TickDir', 'out');  % Set tick marks to face out
% %% plot the first PC
% plot(pc_roc(:,1));
% box off
% set(gca, 'TickDir', 'out');  % Set tick marks to face out

%
% Define a colormap for the subjects
subject_colors = [
    0 0 0;  % Dusty Pink 
    0.15 0.15 0.15;  % Milkshake Purple
    0.3 .3 .3;  % Creamsicle Orange
    0.45 0.45 0.45; % Cream White
    0.6 0.6 0.6;
    0.75 0.75 0.75;
];  % Normalize to [0, 1]


standard_font_size = 32;numbers_font_size = 20;

% Plot the ROC values
subplot('Position', [0.05, 0.1, 0.55, 0.8]);
imagesc(sorted_roc_tc);
hold on
colormap bone;
title('ROC values', 'FontSize', standard_font_size);
xlabel('Time (bins)', 'FontSize', standard_font_size);
ylabel('Electrodes', 'FontSize', standard_font_size);
set(gca, 'FontSize', numbers_font_size);  % Set font size for axis numbers
set(gca, 'TickDir', 'out');  % Set tick marks to face out

if ~isnan(pc_roc)
    clustlines = nan(3, nclust);
    for i = 1:nclust
        temp = sum(clustlabel <= i) + 0.5;
        clustlines(1:2, i) = [temp; temp];
    end
    clustlines = clustlines(:);

    %imagesc(x, y, neurontc(plotorder, :));
    %axis(gca, 'tight', 'ij');
    xlines = cat(1, repmat(xlim', [1, nclust]), nan(1, nclust));
    xlines = xlines(:);
    plot(xlines, clustlines, 'r');
else
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
subplot('Position', [0.62, 0.1, 0.05, 0.8]);
imagesc(subject_color_block);
title('Subjects', 'FontSize', standard_font_size);
yticks([]);  % Remove Y-axis tick marks
set(gca, 'FontSize', numbers_font_size);  % Set font size for axis numbers
set(gca, 'TickDir', 'out');  % Set tick marks to face out

% Create a color block for the coordinates using the thermal colormap
jet_colormap = spring(256);  % Thermal colormap

% Normalize the coordinates to the range [1, 256] for indexing into the colormap
norm_coordinates = round(rescale(sorted_coordinates, 1, 256));

% Create the color block for the coordinates
coordinate_color_block = jet_colormap(norm_coordinates, :);

% Plot the coordinate color block
subplot('Position', [0.69, 0.1, 0.05, 0.8]);
imagesc(reshape(coordinate_color_block, length(sorted_coordinates), 1, 3));
title('Coordinates', 'FontSize', standard_font_size);
yticks(1:10:length(sorted_coordinates));
yticklabels(arrayfun(@num2str, sorted_coordinates(1:10:end), 'UniformOutput', false));
set(gca, 'FontSize', numbers_font_size);  % Set font size for axis numbers
set(gca, 'TickDir', 'out');  % Set tick marks to face out

% Ensure the entire figure fits into the view
set(gcf, 'Position', [100, 100, 1200, 600]);  % Adjust the figure size as needed

% Adjust the position of the coordinates color block subplot
pos = get(gca, 'Position');
pos(1) = 0.75; % Move to the right
pos(3) = 0.05; % Adjust the width
set(gca, 'Position', pos);

% Add the legend for subjects to the right of the coordinates color block
subplot('Position', [0.84, 0.1, 0.05, 0.8]);
axis off;
hold on;
for i = 1:length(folders)
    plot(NaN, NaN, 's', 'MarkerFaceColor', subject_colors(i, :), 'MarkerEdgeColor', 'none');
end
legend(folders, 'Location', 'eastoutside', 'FontSize', numbers_font_size);
hold off;


%% look at group PSTH

figure(2)

clust_indices = find(clustlabel==2);

t = 2;num_conditions = 2;

signal_bounds = [3000 7000];
baseline_bounds = [4800 5000];
window = signal_bounds(2) - signal_bounds(1) + 1;
smoothing_kernel = 1;
x2 = round(signal_bounds(1):smoothing_kernel:signal_bounds(2));
colorflag = 0;

for j = 1:num_conditions

    for_psth = zeros(length(clust_indices), window);

    for i = 1:length(clust_indices)

        n = flagged_subject_ID(clust_indices(i));
        e = flagged_electrode_ID(clust_indices(i));

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


[h, p, ci, stats]=ttest([sum(for_anova{2}(:,1700:2000),2)])
%% plot a heatmap of probability of being in particular cluster and extract data for chi-square test

which_cluster =  3;
% Parameters
window_size = 10;  % Size of the sliding window in mm
overlap = 5;    % Overlap between windows in mm

% Initialize lists for x and y coordinates
all_med_coords = [];
all_trans_coords = [];
cluster_labels = [];

% Collect all coordinates and cluster labels
for i = 1:length(clustlabel)
    n = flagged_subject_ID(i);
    e = flagged_electrode_ID(i);
    med_coor = subject(n).electrode(e).med_coor;
    trans_coor = subject(n).electrode(e).trans_coor;
    
    all_med_coords = [all_med_coords; med_coor];
    all_trans_coords = [all_trans_coords; trans_coor];
    cluster_labels = [cluster_labels; clustlabel(i)];
    %cluster_labels = [cluster_labels; tag(i)];
end

% Determine the min and max coordinates for the grid
x_min = floor(min(all_med_coords))-window_size;
x_max = ceil(max(all_med_coords))+window_size;
y_min = floor(min(all_trans_coords))-window_size;
y_max = ceil(max(all_trans_coords))+window_size;

% Generate sliding window centers with overlap
x_centers = x_min:overlap:x_max;
y_centers = y_min:overlap:y_max;

% Initialize matrix to store the number of Cluster 1 and total electrodes in each bin
cluster1_electrodes_in_bins = zeros(length(y_centers), length(x_centers));
total_electrodes_in_bins = zeros(length(y_centers), length(x_centers));

% Loop over sliding window centers to calculate the number of Cluster 1 and total electrodes
for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        % Define the window range
        x_start = x_centers(x_idx) - window_size/2;
        x_end = x_centers(x_idx) + window_size/2;
        y_start = y_centers(y_idx) - window_size/2;
        y_end = y_centers(y_idx) + window_size/2;
        
        % Find electrodes within the current window
        in_window = all_med_coords >= x_start & all_med_coords < x_end & ...
                    all_trans_coords >= y_start & all_trans_coords < y_end;
        
        % Count total electrodes and Cluster 1 electrodes
        total_electrodes = sum(in_window);
        cluster1_electrodes = sum(in_window & cluster_labels == which_cluster);
        
        % Store the counts
        total_electrodes_in_bins(y_idx, x_idx) = total_electrodes;
        cluster1_electrodes_in_bins(y_idx, x_idx) = cluster1_electrodes;
    end
end

% Calculate the overall expected proportion of Cluster 1 electrodes (null hypothesis)
total_cluster1 = sum(cluster_labels == which_cluster);
total_electrodes_overall = length(cluster_labels);
expected_proportion_cluster1 = total_cluster1 / total_electrodes_overall;

% Calculate the expected number of Cluster 1 electrodes in each bin under the null hypothesis
expected_cluster1_in_bins = total_electrodes_in_bins * expected_proportion_cluster1;

% Reshape the matrices into vectors for chi-square test
observed_cluster1 = cluster1_electrodes_in_bins(:);
expected_cluster1 = expected_cluster1_in_bins(:);

% Remove NaN or zero entries (bins with no electrodes)
valid_bins = total_electrodes_in_bins(:) > 0;
observed_cluster1 = observed_cluster1(valid_bins);
expected_cluster1 = expected_cluster1(valid_bins);

% Run the chi-square square test
% Calculate chi-square statistic manually
chi2_stat = sum((observed_cluster1 - expected_cluster1).^2 ./ expected_cluster1);

% Degrees of freedom: number of bins - 1
df = length(observed_cluster1) - 1;

% Compute p-value
p_value = 1 - chi2cdf(chi2_stat, df);

% Display results
fprintf('Chi-Square Statistic: %.2f\n', chi2_stat);
fprintf('Degrees of Freedom: %d\n', df);
fprintf('p-value: %.4f\n', p_value);


% Display the result of the chi-square test
if p_value>0.05
    disp('The observed distribution is not significantly different from the expected distribution.');
else
    disp('The observed distribution is significantly different from the expected distribution.');
end
disp(['p-value: ', num2str(p_value)]);

% Custom colormap: hot pink -> gray -> blue
cmap_size = 256;  % Size of the colormap
hot_pink = [0 .8 0];  % Hot pink
gray = [0.5, 0.5, 0.5];  % Gray (for null probability)
blue = [0.5, 0.3, 0.8];        % Blue (for lowest probabilities)

% Create the colormap and ensure values are within [0, 1]
custom_cmap = zeros(cmap_size, 3);
for i = 1:cmap_size
    if i <= cmap_size / 2
        % Transition from blue to gray
        ratio = (i - 1) / (cmap_size / 2 - 1);
        custom_cmap(i, :) = blue * (1 - ratio) + gray * ratio;
    else
        % Transition from gray to hot pink
        ratio = (i - cmap_size / 2) / (cmap_size / 2 - 1);
        custom_cmap(i, :) = gray * (1 - ratio) + hot_pink * ratio;
    end
end

% Ensure all values in the colormap are within the range [0, 1]
custom_cmap = max(min(custom_cmap, 1), 0);

% Calculate the null probability
null_probability = total_cluster1 / total_electrodes_overall;

% Set bins with no electrodes to NaN for white background
plot_data = cluster1_electrodes_in_bins ./ total_electrodes_in_bins;
plot_data(total_electrodes_in_bins == 0) = NaN;  % Set empty bins to NaN

% Plot the heatmap
figure;
h = imagesc(x_centers, y_centers, plot_data);  % Plot the data
colormap(custom_cmap);  % Apply the custom colormap
colorbar;
clim([0, 1]);  % Set the color limits to 0 to 1 (probability range)

% Set transparency based on the number of electrodes in each bin
alpha_data = total_electrodes_in_bins./max(max(total_electrodes_in_bins));
%alpha_data = (total_electrodes_in_bins>0);%./max(max(total_electrodes_in_bins));
set(h, 'AlphaData', alpha_data);

% Set the null probability to correspond to gray in the colormap
set(gca, 'CLim', [0, 1]);  % Color limits between 0 and 1 (for probability)
%caxis([0 1]);  % Set color axis with gray for null probability


xlim([-20 40])
ylim([-20 25])

% Ensure y-axis is not flipped
set(gca, 'YDir', 'normal'); 

% Add labels and title
title('Observed Cluster 3 Electrode Distribution', 'FontSize', 14);
xlabel('Medial/Lateral (mm)', 'Fontsize', 14);
ylabel('Rostral/Caudal (mm)', 'Fontsize', 14);
set(gca, 'TickDir', 'out');
box off;
%% ask if cluster 1 and cluster 2 are different (must run the code above with cluster 1 for this to work!)
which_cluster = 2;
% Parameters
window_size = 10;  % Size of the sliding window in mm
overlap = 5;    % Overlap between windows in mm

% Initialize lists for x and y coordinates
all_med_coords = [];
all_trans_coords = [];
cluster_labels = [];

% Collect all coordinates and cluster labels
for i = 1:length(clustlabel)
    n = flagged_subject_ID(i);
    e = flagged_electrode_ID(i);
    med_coor = subject(n).electrode(e).med_coor;
    trans_coor = subject(n).electrode(e).trans_coor;
    
    all_med_coords = [all_med_coords; med_coor];
    all_trans_coords = [all_trans_coords; trans_coor];
    cluster_labels = [cluster_labels; clustlabel(i)];
    %cluster_labels = [cluster_labels; tag(i)];
end


% Initialize matrix to store the number of Cluster 1 and total electrodes in each bin
cluster2_electrodes_in_bins = zeros(length(y_centers), length(x_centers));
total_electrodes_in_bins = zeros(length(y_centers), length(x_centers));

% Loop over sliding window centers to calculate the number of Cluster 1 and total electrodes
for x_idx = 1:length(x_centers)
    for y_idx = 1:length(y_centers)
        % Define the window range
        x_start = x_centers(x_idx) - window_size/2;
        x_end = x_centers(x_idx) + window_size/2;
        y_start = y_centers(y_idx) - window_size/2;
        y_end = y_centers(y_idx) + window_size/2;
        
        % Find electrodes within the current window
        in_window = all_med_coords >= x_start & all_med_coords < x_end & ...
                    all_trans_coords >= y_start & all_trans_coords < y_end;
        
        % Count total electrodes and Cluster 1 electrodes
        total_electrodes = sum(in_window);
        cluster2_electrodes = sum(in_window & cluster_labels == 2);
        
        % Store the counts
        total_electrodes_in_bins(y_idx, x_idx) = total_electrodes;
        cluster2_electrodes_in_bins(y_idx, x_idx) = cluster2_electrodes;
    end
end


% Reshape the matrices into vectors for chi-square test
observed_cluster1 = cluster1_electrodes_in_bins(:)+0.1;
expected_cluster2 = cluster2_electrodes_in_bins(:)+0.1;

% Remove NaN or zero entries (bins with no electrodes)
valid_bins = (expected_cluster2) > 0;
observed_cluster1 = observed_cluster1(valid_bins);
expected_cluster2 = expected_cluster2(valid_bins);

% Normalize the expected distribution to have the same total as the observed distribution
expected_cluster2_normalized = expected_cluster2 * (sum(observed_cluster1) / sum(expected_cluster2));

% Run the chi-square test
chi2_stat = sum((observed_cluster1 - expected_cluster2_normalized).^2 ./ expected_cluster2_normalized);

% Degrees of freedom: number of bins - 1
df = length(observed_cluster1) - 1;

% Compute p-value
p_value = 1 - chi2cdf(chi2_stat, df);

% Display results
fprintf('Chi-Square Statistic: %.2f\n', chi2_stat);
fprintf('Degrees of Freedom: %d\n', df);
fprintf('p-value: %.4f\n', p_value);

% Display the result of the chi-square test
if p_value > 0.05
    disp('The distributions of Cluster 1 and Cluster 2 are not significantly different.');
else
    disp('The distributions of Cluster 1 and Cluster 2 are significantly different.');
end

