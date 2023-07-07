function [] = calculate_z_score_w_consistent_baselines()
% Previously I calculate z score using baselines that are different for
% each stimulus. Now I will pool all calculated baselines to have more
% comparable z scores.

    %% Load signal events
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
            if contains(child_name, 'Entire')
                continue
            end
            load(horzcat(file, name, '/', child_name, '/dFF_data.mat'));
            all_baselines = horzcat(all_baselines, reshape(baselines, [1, (size(baselines, 1)*size(baselines, 2))]));
        end
        %Calculate std_avg for all data
        std_avg = std(all_baselines);
        avg = mean(all_baselines);

        for child_ind = 3:size(children, 1)
            child = children(child_ind);
            child_name = child.name;
            if contains(child_name, 'Entire')
                continue
            end
            load(horzcat(file, name, '/', child_name, '/dFF_data.mat'));
            diffs = (US_trials_percent - avg);
            z_US_trials = (diffs ./ std_avg);
            save(horzcat(file, name, '/', child_name, '/consistent_z_scores.mat'), 'z_US_trials');

        end


    end
end
