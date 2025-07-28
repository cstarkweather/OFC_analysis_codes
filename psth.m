function psth(trajectories_x,color)

    % Suppose trajectories_x is size [m x n].
    %  m = number of trials
    %  n = number of time points

    [m, n] = size(trajectories_x);

    % If you have a time vector, use it here. Otherwise, define one:
    time = 1:n;

    % 1) Compute mean and std across trials (i.e., across rows)
    meanTrace = mean(trajectories_x, 1);  % size [1 x n]
    stdTrace  = std( trajectories_x, 0, 1);  % size [1 x n], 0 means sample-based stdev

    % 2) Build upper and lower bands for plotting
    upperBand = meanTrace + stdTrace/sqrt(m) *2;
    lowerBand = meanTrace - stdTrace/sqrt(m) *2;


    % (A) Plot the shaded region for +/- 1 std
    %  We'll use 'fill' to create a patch:
    fill([time, fliplr(time)], ...
         [upperBand, fliplr(lowerBand)], ...
         color, ...            % color = black
         'FaceAlpha', 0.2, ...  % transparency
         'EdgeColor','none');   % no border line
hold on
    % (B) Plot the mean trace on top
    plot(time, meanTrace, '-', 'Color',color,'LineWidth', 2);

    xlabel('Time (samples)');
    ylabel('Amplitude');
    title('Mean \pm 1SD Across Trials');
    box off;
    set(gca, 'TickDir', 'out');
    hold on
    
end
