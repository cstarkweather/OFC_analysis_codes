function binned_data = bin_data_columns(bin_width, trialIDs, signal)
    % Convert bin width to number of columns
    n_columns = bin_width;

    % Determine the number of bins and trials
    n_bins = ceil(size(signal, 2) / n_columns);

    % Calculate the number of columns to pad
    n_pad_columns = n_bins * n_columns - size(signal, 2);

    % Pad the signal matrix to have a multiple of n_columns
    padded_signal = [signal, NaN(size(signal, 1), n_pad_columns)];

    % Reshape the padded signal to group columns into bins
    reshaped_signal = reshape(padded_signal(:, 1:n_bins * n_columns), size(signal, 1), n_columns, n_bins);

    % Compute the mean for each bin along the second dimension (column axis)
    binned_data = squeeze(nanmean(reshaped_signal(trialIDs, :, :), 2));
end