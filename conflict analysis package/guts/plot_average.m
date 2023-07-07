function [] = plot_average(mode, brain_area, uniform_trials, seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, laser_duration)
    %Plot average dF/F signal in %. Takes an N by L matrix
    %TODO: take out NANs
    
    CS_duration = 1;
    
    US_at_o = 0;
    for row = 1:size(uniform_trials, 1)
        for col = 1:size(uniform_trials, 2)
            if isnan(uniform_trials(row,col))
                uniform_trials(row,col) = 0;
            end
        end
    end
    
    uniform_sum_signal = sum(uniform_trials, 1);
    
    uniform_average_signal = uniform_sum_signal/size(uniform_trials, 1);
    
    smallest_size = size(uniform_trials, 2);
    why = figure;
    hold on
    x_axis = ((1:smallest_size)/smallest_size) * seconds_per_trial; %Divide into 6 seconds, half less than 0
    x_axis = x_axis - x_axis(US_target_frame);
    
    time_of_US = 0;
    time_of_CS = 0;
    bigness = size(uniform_average_signal);
    
    if US_at_o
        x_axis = x_axis - time_of_US;
        time_of_US = 0;
    end
    
    
    if and(sound_only, string(mode) == '110')
        plot(x_axis, uniform_average_signal);
        yLimits = get(gca,'YLim'); 
        Water_CS = plot([time_of_CS, time_of_CS], yLimits, 'color', 'r', 'DisplayName','CS', 'MarkerSize', 20);
        legend([Water_CS], {'Water CS'});
    elseif or(strcmp(alignment, 'laser'), strcmp(alignment, 'events'))
        plot(x_axis, uniform_average_signal);
        yLimits = get(gca,'YLim');
        hold off
        close all
        figure
        hold on
        x1 = [time_of_US, time_of_US + laser_duration];
        y1 = [yLimits(1), yLimits(1)];%[numel(RawEvent),numel(RawEvent)];
        Laser=area(x1, y1,'LineStyle','None');
        Laser(1).FaceColor = [0.976 0.894 0.690];
        x1 = [time_of_US, time_of_US + laser_duration];
        y1 = [yLimits(2), yLimits(2)];%[numel(RawEvent),numel(RawEvent)];
        Laser=area(x1, y1,'LineStyle','None');
        Laser(1).FaceColor = [0.976 0.894 0.690];
        plot(x_axis, uniform_average_signal);
        set(gca, 'ylim',yLimits,'FontSize', 10);
        legend(Laser, {alignment})
    elseif strcmp(alignment, 'avoid')
        plot(x_axis, uniform_average_signal);
        yLimits = get(gca,'YLim');
        Avoid = plot([time_of_US, time_of_US], yLimits, 'color', 'c', 'DisplayName','Avoidance Initiation', 'MarkerSize', 20);
        legend(Avoid, {'Avoidance Initiation'})
    elseif strcmp(alignment, 'sound')
        plot(x_axis, uniform_average_signal);
        yLimits = get(gca,'YLim');
        so = plot([time_of_CS, time_of_CS], yLimits, 'color', 'c', 'DisplayName','CS', 'MarkerSize', 20);
        legend(so, {'CS'})
    elseif strcmp(alignment, 'baseline run')
        plot(x_axis, uniform_average_signal);
        yLimits = get(gca,'YLim');
        set(gca, 'xlim',[0 4]);
        %yLimits = get(gca,'YLim');
        Avoid = plot([time_of_US, time_of_US], yLimits, 'color', 'c', 'DisplayName','Running Initiation', 'MarkerSize', 20);
        legend(Avoid, {'Running Initiation'})
    elseif strcmp(alignment, 'escape')
        plot(x_axis, uniform_average_signal);
        yLimits = get(gca,'YLim');
        set(gca, 'xlim',[4 10]);
        %yLimits = get(gca,'YLim');
        Avoid = plot([time_of_US, time_of_US], yLimits, 'color', 'c', 'DisplayName','Running Initiation', 'MarkerSize', 20);
        legend(Avoid, {'Escape'})
    elseif and(sound_only, mode == '101')
        plot(x_axis, uniform_average_signal);
        yLimits = get(gca,'YLim');
        %yLimits = get(gca,'YLim');
        x1 = [time_of_CS,time_of_CS + CS_duration];
        y1 = yLimits;%[numel(RawEvent),numel(RawEvent)];
        %y1 = [numel(RawEvent),numel(RawEvent)];
        h=brain_area(x1,y1,'LineStyle','None');
        h(1).FaceColor = [0.976 0.894 0.690];
        %Air_CS = plot([time_of_CS, time_of_CS], yLimits, 'color', 'r', 'DisplayName','CS', 'MarkerSize', 20);
        legend([Air_CS], {'Air CS'});
    elseif mode == '110'
        hold off
        get_y = plot(x_axis, uniform_average_signal);
        yLimits = get(gca,'YLim');
        %hold off
        close(why);
        figure
        hold on
        
        %Plot the US 
        x1 = [time_of_US, time_of_US + .2];
        y1 = [yLimits(1), yLimits(1)];%[numel(RawEvent),numel(RawEvent)];
        US_period=area(x1, y1,'LineStyle','None');
        US_period(1).FaceColor = [153/255 1 1];
        x1 = [time_of_US, time_of_US + .2];
        y1 = [yLimits(2), yLimits(2)];%[numel(RawEvent),numel(RawEvent)];
        US_period=area(x1, y1,'LineStyle','None');
        US_period(1).FaceColor = [153/255 1 1];
        
        %Plot the CS 
        x1 = [time_of_CS, time_of_CS + laser_duration];
        y1 = [yLimits(1), yLimits(1)];%[numel(RawEvent),numel(RawEvent)];
        CS_period=area(x1, y1,'LineStyle','None');
        CS_period(1).FaceColor = [0.976 0.894 0.690];
        x1 = [time_of_CS, time_of_CS + laser_duration];
        y1 = [yLimits(2), yLimits(2)];%[numel(RawEvent),numel(RawEvent)];
        CS_period=area(x1, y1,'LineStyle','None');
        CS_period(1).FaceColor = [0.976 0.894 0.690];
        
        
        plot(x_axis, uniform_average_signal);
        set(gca, 'ylim',yLimits,'FontSize', 10, 'xlim', [-10, 10]);
        legend([US_period, CS_period], {'Water', 'CS'});
    elseif mode == '100'
        plot(x_axis, uniform_average_signal);
        %yLimits = get(gca,'YLim');
        yLimits = get(gca,'YLim');
        CS = plot([time_of_CS time_of_CS], yLimits, 'color', 'r', 'DisplayName','CS', 'MarkerSize', 20);
        legend(CS, {'CS'});
    elseif mode == '010'
        plot(x_axis, uniform_average_signal);
        %yLimits = get(gca,'YLim');
        yLimits = get(gca,'YLim');
        Water = plot([time_of_US, time_of_US], yLimits, 'color', 'c', 'DisplayName','US', 'MarkerSize', 20);
        legend(Water, {'Water'})
    elseif mode == '001'
        plot(x_axis, uniform_average_signal);
        %yLimits = get(gca,'YLim');
        yLimits = get(gca,'YLim');
        Air = plot([time_of_US, time_of_US], yLimits, 'color', 'c', 'DisplayName','US', 'MarkerSize', 20);
        legend(Air, {'Air'})
    elseif mode == '101'
        %{
        plot(x_axis, uniform_average_signal);
        %yLimits = get(gca,'YLim');
        yLimits = get(gca,'YLim');
        Air = plot([time_of_US, time_of_US], yLimits, 'color', 'c', 'DisplayName','US', 'MarkerSize', 20);
        CS = plot([time_of_CS, time_of_CS], yLimits, 'color', 'r', 'DisplayName','CS', 'MarkerSize', 20);
        legend([Air, CS], {'Air', 'CS'})
        %}
        
        
        
        get_y = plot(x_axis, uniform_average_signal);
        yLimits = get(gca,'YLim');
        hold off
        close(get_y);
        figure
        hold on
        
        %Plot the US 
        x1 = [time_of_US, time_of_US + .2];
        y1 = [yLimits(1), yLimits(1)];%[numel(RawEvent),numel(RawEvent)];
        US_period=area(x1, y1,'LineStyle','None');
        US_period(1).FaceColor = [252/255 189/255 189/255];
        x1 = [time_of_US, time_of_US + .2];
        y1 = [yLimits(2), yLimits(2)];%[numel(RawEvent),numel(RawEvent)];
        US_period=area(x1, y1,'LineStyle','None');
        US_period(1).FaceColor = [252/255 189/255 189/255];
        
        %Plot the CS 
        x1 = [time_of_CS, time_of_CS + laser_duration];
        y1 = [yLimits(1), yLimits(1)];%[numel(RawEvent),numel(RawEvent)];
        CS_period=area(x1, y1,'LineStyle','None');
        CS_period(1).FaceColor = [0.976 0.894 0.690];
        x1 = [time_of_CS, time_of_CS + laser_duration];
        y1 = [yLimits(2), yLimits(2)];%[numel(RawEvent),numel(RawEvent)];
        CS_period=area(x1, y1,'LineStyle','None');
        CS_period(1).FaceColor = [0.976 0.894 0.690];
        
        
        plot(x_axis, uniform_average_signal);
        set(gca, 'ylim',yLimits,'FontSize', 10);
        legend([US_period, CS_period], {'Air', 'CS'});
        
        
        
        
    else
        disp('OH NO')
    end

    title(horzcat('DA release in ', brain_area, ', average % dF/F'), 'FontSize', 20);
    xlabel('Seconds');
    ylabel('dF/F');
    %xlim([0 10]);
    %ylim([-1 2]);
    hold off

end

