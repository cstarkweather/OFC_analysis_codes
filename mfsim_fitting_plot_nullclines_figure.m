function [xx, yy, tt] = mfsim_fitting_plot_nullclines_figure( f1, f2, pausetime, doplot,alphaReward,alphaPunish,noise,lambda,w,alpha,valence,offset,visualize_trajectory)

    xgrd = -5000:0.1:5000;       % i/o function: x
    tau  = 40;    % time constant of dynamics

    t    = 1:6000;     % note: time step dt = 1 msec

    %--------------external inputs
    %--------------------------------------
    % Example: Single-phase approach/avoid
    %--------------------------------------

    % Clear the old input definitions
    Ex = zeros(1, length(t));
    Ey = zeros(1, length(t));

    % Suppose the subject sees bombs & treasure from t=1ms to t=6000ms
    tStart = 1;
    tEnd   = 6000;


    baselineApproach = 0;   % baseline drive to approach
    baselineAvoid    = 0;   % baseline drive to avoid

    % For the entire trial, just set a baseline
    Ex(:) = 0;
    Ey(:) = 0;

    % During tStart to tEnd, add the relevant magnitude
    %   (example: approach reward minus punishment, etc.)
    Ex(tStart:tEnd) = Ex(tStart:tEnd) + alphaReward*f1  -alphaPunish*f2;
    Ey(tStart:tEnd) = Ey(tStart:tEnd) ;




    %--------------initialize plots
    pathHandle = [];  % Will be the thin black line showing the path
    posxy      = [];  % Will be the red dot
    ttHandle   = [];  % Text handle

    if (doplot)

        %figure(1); clf;
        pos = get(gcf,'Position');
        set(gcf,'Position', [pos(1) pos(2) 360 600] );

        %------------------ Top subplot ------------------
        %subplot('Position', [0.15 0.5 0.7 0.4]);
        %figure(1)
        hold on;
        % The nullclines you’re already plotting
        nc1 = plot( xgrd, xgrd, 'k' );
        nc2 = plot( xgrd, xgrd, 'k' );

        % (Below lines used to have EraseMode='xor', but we remove that for modern usage)
        set(nc1, 'Color', 'k','Linewidth',2);
        set(nc2, 'Color', 'k','Linewidth',2);

        % 1) Create a thin black line for the path (initially empty)
        pathHandle = plot(NaN, NaN, 'k-', 'LineWidth', 1);

        % 2) Create the red dot for the current position
        posxy = plot(0, 0, 'r.', 'MarkerSize', 20);

        xlabel('nS (plus neuron)');
        ylabel('nS (minus neuron)');
        axis( [-0.5 4 -0.5 4]);
        set(gca, 'TickDir', 'out');

        % A text label for time
        ttHandle = text(3, 4.5, 't=0');

        %------------------ Bottom subplot ------------------
        %subplot('Position', [0.15 0.15 0.7 0.2]);
        %traj = plot(0, 0, 'r');
        %xlabel('time (msec)');
        %ylabel('nS (plus neuron)');
        % axis( [0 par.T 0 5] );
        %set(gca, 'TickDir', 'out');

    end

    %--------------simulation
    x = 0; 
    y = 25;
    xx = []; 
    yy = [];

    % These will store the position history for the black line
    xPath = [];
    yPath = [];


    for k = 1:length(t)

        % ODE integration using Euler's method
        dx = -x + iofunc( -w*y + lambda*Ex(k),alpha,-valence,-offset) + noise*randn;
        dy = -y + iofunc( -w*x + lambda*Ey(k),alpha,-valence,-offset) + noise*randn;

        x = x + dx/tau;   % dt=1, so dx*dt = dx
        y = y + dy/tau;

        % Every 100th iteration, update the plot

        if doplot
            samplefactor = 2; % downsample the trajectory by a factor of 2 if plotting
        else
            samplefactor = 1;
        end

        if ( ~mod(k,samplefactor) )
            xx = [xx, w*x]; 
            yy = [yy, w*y];

            %================= Update the visuals =================
            if (doplot)
                % 1) Update the nullclines (scaled by 1000, as in your original code)
                set(nc1, 'XData', xgrd, ...
                         'YData', w*iofunc( -xgrd + lambda*Ey(k),alpha,-valence,-offset));
                set(nc2, 'XData', w*iofunc( -xgrd + lambda*Ex(k),alpha,-valence,-offset), ...
                         'YData', xgrd );

                % 2) Append this (x,y) to the path and update
                %xPath = [xPath, 1000*w*x];
                %yPath = [yPath, 1000*w*y];
                %set(pathHandle, 'XData', xPath, 'YData', yPath);

                % 3) Move the red dot to the current position
                %set(posxy, 'XData', 1000*w*x, 'YData', 1000*w*y);


                % 5) Update the bottom subplot (time vs plus neuron)
                %set(traj, 'XData', (1:k/100)*10, 'YData', 1000*xx);

                %drawnow;
                %pause( pausetime ); % slows down the dynamics for visualization
            end
        end
    end

    % Return time indices for the downsampled path
    tt = t(1:samplefactor:end);
if doplot
% Pick a single time in [tStart, tEnd]
kSnapshot = round((tStart + tEnd)/2);

% Build a grid in the same scale as your top subplot
nGrid = 9;
xVals = linspace(-0.5,4, nGrid);  % or [-350*2, 350*2]
yVals = linspace(-0.5,4, nGrid);
[Xg, Yg] = meshgrid(xVals, yVals);
DX = zeros(size(Xg));
DY = zeros(size(Xg));

% Evaluate the input terms at that snapshot
ExVal = Ex(kSnapshot); 
EyVal = Ey(kSnapshot);

for i = 1:numel(Xg)
    % Convert from plotted coords to "model" coords
    x_model = Xg(i) / w;
    y_model = Yg(i) / w;

    dx = -x_model + iofunc(-w*y_model + lambda*ExVal,alpha,-valence,-offset);
    dy = -y_model + iofunc(-w*x_model + lambda*EyVal,alpha,-valence,-offset);

    % Scale them back to plot coords
    DX(i) = dx*w*10;
    DY(i) = dy*w*10;
end

% Now plot these arrows on the top subplot
% figure(1);
% subplot('Position', [0.15 0.5 0.7 0.4]);
hold on;
quiver(Xg, Yg, DX, DY, 'Color','b', 'AutoScale','on','AutoScaleFactor',1);
else
    
end

box off


    function y = iofunc(x, alpha,valence,offset)
    %alpha is an overall gain/scale,
    %beta controls smoothness around x=0 (larger beta -> steeper).
    y = alpha * tanh(x-valence) + offset;
    %y = valence * alpha * log(1 + exp(valence * beta * (x + offset))) / beta ;
end
end
   %y = alpha * max(0, x);
   %y = alpha * tanh(x);
