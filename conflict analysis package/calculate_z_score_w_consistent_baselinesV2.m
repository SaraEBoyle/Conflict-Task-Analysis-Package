function [] = calculate_z_score_w_consistent_baselinesV2()
% Previously I calculate z score using baselines that are different for
% each stimulus. Now I will pool all calculated baselines to have more
% comparable z scores.

    %% Load signal events
    left_trials = [];
    right_trials = [];
    file = horzcat(pwd, '/');
    tem = dir(horzcat(file, '*', 'signal')); 
    for signal_folder_ind = 1:size(tem, 1)
        signal_folder = tem(signal_folder_ind, :);
        name = signal_folder.name;
        children = dir(horzcat(file, name, '/')); 
        all_baselines = [];
        for child_ind = 3:size(children, 1)
            child = children(child_ind);
            child_name = child.name;
            if or(contains(child_name, 'Entire'), contains(child_name, 'pearson'))
                continue
            elseif contains(child_name, 'full_session')
                continue
            elseif contains(child_name, 'raw_data')
                continue
            end
            load(horzcat(file, name, '/', child_name, '/dFF_data.mat'));
            all_baselines = horzcat(all_baselines, reshape(baselines, [1, (size(baselines, 1)*size(baselines, 2))]));
        end
        %Calculate std_avg for all data
        std_avg = std(all_baselines);
        avg = mean(all_baselines);

        for child_ind = 3:size(children, 1)
            US_trials = [];
            child = children(child_ind);
            child_name = child.name;
            if or(contains(child_name, 'Entire'), contains(child_name, 'pearson'))
                continue
            elseif contains(child_name, 'full_session')
                continue
            elseif contains(child_name, 'raw_data')
                continue
            end
            load(horzcat(file, name, '/', child_name, '/dFF_data.mat'));
            if ~isempty(left_trials)
                US_trials_percent = left_trials;
                left_trials = [];
            elseif ~isempty(right_trials)
                US_trials_percent = right_trials;
                right_trials = [];
            elseif ~isempty(US_trials)
                US_trials_percent = US_trials;
                US_trials = [];
            end
            diffs = (US_trials_percent - avg);
            z_US_trials = (diffs ./ std_avg);
            cur_file = horzcat(file, name, '/', child_name);
            if exist(horzcat(cur_file, '/artifacts_removed'), 'dir')
                load(horzcat(cur_file, '/', 'artifacts_removed/deleted_trials.mat'));
                z_US_trials(to_remove, :) = [];
            end
            %mean_z_US_trials = mean(z_US_trials);
            
            if size(z_US_trials, 1) > 1
                mean_z_US_trials = mean(z_US_trials);
            else
                mean_z_US_trials = z_US_trials; 
            end
            lefts = find(x_axis >= 0);
            rights = find(x_axis <= 5);
            zero_thru_five = intersect(lefts, rights);
            zero_thru_five_axis = x_axis(zero_thru_five);
            zero_thru_five_mean_z = mean_z_US_trials(zero_thru_five);
            zero_thru_five_AUC = trapz(zero_thru_five_axis, zero_thru_five_mean_z);
            
            save(horzcat(file, name, '/', child_name, '/consistent_z_scores.mat'), 'z_US_trials', 'mean_z_US_trials', 'zero_thru_five_AUC');
        end
    end
    
    %% Now for control
    tem = dir(horzcat(file, '*', 'control')); 
    for signal_folder_ind = 1:size(tem, 1)
        signal_folder = tem(signal_folder_ind, :);
        name = signal_folder.name;
        children = dir(horzcat(file, name, '/')); 
        all_baselines = [];
        for child_ind = 3:size(children, 1)
            child = children(child_ind);
            child_name = child.name;
            if or(contains(child_name, 'Entire'), contains(child_name, 'pearson'))
                continue
            end
            load(horzcat(file, name, '/', child_name, '/dFF_data.mat'));
            all_baselines = horzcat(all_baselines, reshape(baselines, [1, (size(baselines, 1)*size(baselines, 2))]));
        end
        %Calculate std_avg for all data
        std_avg = std(all_baselines);
        avg = mean(all_baselines);

        for child_ind = 3:size(children, 1)
            US_trials = [];
            child = children(child_ind);
            child_name = child.name;
            if or(contains(child_name, 'Entire'), contains(child_name, 'pearson'))
                continue
            elseif contains(child_name, 'full_session')
                continue
            elseif contains(child_name, 'raw_data')
                continue
            end
            load(horzcat(file, name, '/', child_name, '/dFF_data.mat'));
            if ~isempty(left_trials)
                US_trials_percent = left_trials;
                left_trials = [];
            elseif ~isempty(right_trials)
                US_trials_percent = right_trials;
                right_trials = [];
            elseif ~isempty(US_trials)
                US_trials_percent = US_trials;
                US_trials = [];
            end
            diffs = (US_trials_percent - avg);
            z_US_trials = (diffs ./ std_avg);
            fake_name = strrep(name, 'control', 'signal');
            cur_file = horzcat(file, fake_name, '/', child_name);
            if exist(horzcat(cur_file, '/artifacts_removed'), 'dir')
                load(horzcat(cur_file, '/', 'artifacts_removed/deleted_trials.mat'));
                z_US_trials(to_remove, :) = [];
            end
            if size(z_US_trials, 1) > 1
                mean_z_US_trials = mean(z_US_trials);
            else
                mean_z_US_trials = z_US_trials; 
            end
            lefts = find(x_axis >= 0);
            rights = find(x_axis <= 5);
            zero_thru_five = intersect(lefts, rights);
            zero_thru_five_axis = x_axis(zero_thru_five);
            zero_thru_five_mean_z = mean_z_US_trials(zero_thru_five);
            zero_thru_five_AUC = trapz(zero_thru_five_axis, zero_thru_five_mean_z);
            save(horzcat(file, name, '/', child_name, '/consistent_z_scores.mat'), 'z_US_trials', 'mean_z_US_trials', 'zero_thru_five_AUC');
        end
    end
    
end
