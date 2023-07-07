function [] = plot_heatmap(uniform_trials, area, seconds_per_trial, US_target_frame, CS_target_frame, type, pre_time)
%UNTITLED2 Summary of this function goes here
%% Plot trial by trial responses
figure;

colormap('parula');
uniform_trials = uniform_trials*100; % Change to percent
imagesc(uniform_trials);
%Use this to set the axis yourself
%caxis([-3 7]);
colorbar;
ten = round(size(uniform_trials, 2)/seconds_per_trial);
x_axis = horzcat(1, (1:seconds_per_trial)*ten);
xticks(x_axis);

title(horzcat('DA release in ', area,', % dF/F'), 'FontSize', 20);
xticklabels(strsplit(num2str((0:seconds_per_trial) - pre_time),' '))
xlabel('Seconds');
ylabel('Trials');
hold on
x_US = US_target_frame*ones(1, 2*size(uniform_trials, 1));
x_CS = CS_target_frame*ones(1, 2*size(uniform_trials, 1));
%y_US = (1:(2*size(uniform_trials, 1)))/2;
plot(x_US, (1:(2*size(uniform_trials, 1)))/2, '.', 'color', 'white');
if and(~strcmp(type, 'baseline run'), ~strcmp(type, 'avoid'))
    plot(x_CS, (1:(2*size(uniform_trials, 1)))/2, '.', 'color', 'c');
elseif strcmp(type, 'baseline run')
    %set(gca, 'xlim', [0 4]);
elseif strcmp(type, 'avoid')
    %set(gca, 'xlim', [0 4]);
    
end
 end

