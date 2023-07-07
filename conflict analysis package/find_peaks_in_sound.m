function [] = find_peaks_in_sound()
    %% Finds the peaks in activity per trial pre sound, during sound, and post sound
    pre_time = 20;
    sound_time = 20;
    post_time = 20;
    pre_data = find(x_axis < 0);
    sound_data_right = find(x_axis > 0);
    sound_data_left = find(x_axis < (sound_time - 2));
    sound_data = intersect(sound_data_left, sound_data_right);
    post_data = find(x_axis > sound_time);
    
    %Split the z_score data into pre, sound, and post
    pre_data = z_US_trials(:, pre_data);
    sound_data = z_US_trials(:, sound_data);
    post_data = z_US_trials(:, post_data);
    
    pre_data_peaks = max(pre_data');
    sound_data_peaks = max(sound_data');
    post_data_peaks = max(post_data');
    
    disp(horzcat('Pre peaks mean: ', num2str(mean(pre_data_peaks))));
    disp(horzcat('Sound peaks mean: ', num2str(mean(sound_data_peaks))));
    disp(horzcat('Post peaks mean: ', num2str(mean(post_data_peaks))));
end