function [decision_time, decision, num_states] =  detectStatePreferences(p,time_stable,thresh)

decision = 2;
for i=1:length(p) - time_stable
    snippet = p(i:i+time_stable-1);
    running_pref = sum(snippet>0.5);
    if (running_pref/time_stable < 1-thresh || running_pref/time_stable > thresh) && (i+time_stable) > 600
        decision_time = i + time_stable;
        decision = running_pref/time_stable > thresh;
        full_run = p(1:(i+time_stable-1));
        %full_run = smooth(full_run,20);
        above_05_curr = full_run(1:end-1) > 0.5;   % whether the current sample is above 0.5
        above_05_next = full_run(2:end)   > 0.5;   % whether the next sample is above 0.5

        crossings_in_window = ( above_05_curr & ~above_05_next ) ...
            | ( ~above_05_curr &  above_05_next );

        num_states = sum(crossings_in_window)+1;

        break
    else
    end
end

if decision ==2
    decision_time = length(p);
    decision = p(end)>0.5;
    full_run = p;
    %full_run = smooth(full_run,20);
    above_05_curr = full_run(1:end-1) > 0.5;   % whether the current sample is above 0.5
    above_05_next = full_run(2:end)   > 0.5;   % whether the next sample is above 0.5

    crossings_in_window = ( above_05_curr & ~above_05_next ) ...
        | ( ~above_05_curr &  above_05_next );

    num_states = sum(crossings_in_window)+1;
else
end

end