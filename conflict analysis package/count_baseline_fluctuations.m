load('event_times.mat')
load('event_times.mat', 'shock_received_times')
load('raw_data.mat');

%Names of fibers you're checking. Order matters here
fiber_names = {'TS', 'CeL'};

%List fibers you want to check in the analysis
signal_fiber = 2;
control_fiber = 2;
data = raw_signal_dFFs{signal_fiber};
control_data = raw_control_dFFs{control_fiber};
if ~exist('baseline activity analysis', 'dir')
    mkdir('baseline activity analysis');
end
file = horzcat('baseline activity analysis/', fiber_names{signal_fiber}, '/');
if ~exist(file, 'dir')
    mkdir(file);
end

%% Plot control channel and find areas of high artifacts
control_fig = figure;
hold on
plot(control_data(:, 1), control_data(:, 2), 'black');
control_threshold = 5*std(control_data(:, 2));
plot([0, 2500], [control_threshold, control_threshold], 'c');
title('control');

control = control_data(:, 2);
time = control_data(:, 1);
above = 0;
artifacts = [];
end_artifacts = [];
for x = 1:size(control, 1)
    if or(control(x) > control_threshold, control(x) < (0 - control_threshold))
        %if over threshold
        if above == 0
            %if the first point above threshold
            point = time(x);
            artifacts = horzcat(artifacts, x);
            above = 1;
        else
            %if not
            continue
        end
    else
        %if under threshold
        if above == 0
            %if previous point under threshold
            continue
        else
            %if previous point above threshold
            above = 0; %Reset
            end_artifacts = horzcat(end_artifacts, x);
        end
        
    end
end
to_delete = [];
for b = 1:size(artifacts, 2)
    ind = artifacts(b) - 5;
    if b > size(end_artifacts, 2)
        to_delete = horzcat(to_delete, ind);
        break
    end
    end_ind = end_artifacts(b) + 5;
    to_delete = horzcat(to_delete, ind:end_ind);
end
to_delete = unique(to_delete);
control_data(to_delete, :) = [];
plot(control_data(:, 1), control_data(:, 2), 'g');
savefig(horzcat(file, '/control_activity_plot.fig'));
hold off
%% delete the high artifact regions from the signal data
new_data = data;
new_data(to_delete, :) = [];
data_fig = figure;
hold on
og_data = plot(data(:, 1), data(:, 2), 'black');
corrected_data = plot(new_data(:, 1), new_data(:, 2), 'g');
top = max(new_data(:, 2));
load('event_times.mat', 'shock_received_times')
p1 = plot([shock_received_times;shock_received_times], [zeros(1, size(shock_received_times, 2));(top * ones(1, size(shock_received_times, 2)))], 'r');
title('signal');

%{
disp('Select interval to analyze. ');
[x_in,y_in] = ginput(2);

%Cut anything funky from the end
left_data = data(find(data(:,1) < x_in(1)), :);
right_data = data(find(data(:,1) > x_in(2)), :);
new_data = vertcat(left_data, right_data);
plot(new_data(:, 1), new_data(:, 2), 'g');
%}
signal = new_data(:, 2);
time = new_data(:, 1);
threshold = 2*std(signal);

p3 = plot([0, 2500], [threshold, threshold], 'c');
legend([og_data, corrected_data, p1(1), p3], 'unedited','removed artifacts', 'shocks', 'high activity threshold');
savefig(horzcat(file, '/signal_activity_plot.fig'));
above = 0;
active = [];
for x = 1:size(signal, 1)
    if signal(x) > threshold
        %if over threshold
        if above == 0
            %if the first point above threshold
            point = time(x);
            for shock = shock_received_times
                if and(point > shock - 1, point < shock + 5)
                    %if the point is during a shock, ignore it
                    above = 1;
                    break
                end
            end
            if above == 1
                continue
            end
            active = horzcat(active, x);
            above = 1;
        else
            %if not
            continue
        end
    else
        %if under threshold
        if above == 0
            %if previous point under threshold
            continue
        else
            %if previous point above threshold
            above = 0; %Reset
            
        end
        
    end
end
active_times = time(active);
save(horzcat(file, '/active_times.mat'), 'active_times');
disp(above);
hold off;
