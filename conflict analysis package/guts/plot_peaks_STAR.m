function [left_peaks, left_mean] = plot_peaks_STAR(corrected_first_trials, first_US_frame, right_window, fps, left_window)
    frames_to_check = right_window*(fps/2);
    left_peaks = [];
    for x = 1:size(corrected_first_trials, 1)
        plot(corrected_first_trials(x,:));
        new_max = max(corrected_first_trials(x,(first_US_frame + left_window * fps/2):(first_US_frame+frames_to_check)));
        hold on
        plot(120, new_max, '*');
        left_peaks = horzcat(left_peaks, new_max);
        hold off
    end
    
    Left = left_peaks;
    left_min = .7;
    left_max = .8;
    left_spread = left_min + (left_max-left_min) .* rand(1, length(left_peaks));
    %x = ones(1, length(Left));
    x = left_spread;
    left_mean = mean(Left);
    
    plot(x, Left, 'bo', 'MarkerSize', 8, 'LineWidth', 3);
    hold on;
    plot(.75, left_mean, 'c*', 'MarkerSize', 20, 'LineWidth', 12);
    xlim([0.5, 1.5]);
    
    ax = gca;
    ax.XTick = [.75, 1.2];
    %grid on;
end