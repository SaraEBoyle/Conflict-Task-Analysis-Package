function [uniform_trials, US_target_frame, true_seconds_per_trial, sec_offsets, true_sec_offsets] = sort_and_align_trials_V7(trials, pre_ITIs, post_ITIs, ...
    seconds_per_trial, US_frame, frame_rate, seconds_to_extend,  do_alignment)
    %Align a set of trial data, given trials size number of trials by
    %frame, previous ITI data, post ITI, baseline frame number, trial by
    %trial correction. Returns aligned all trials, US trials, and no US
    %trials.
    %TODO: sound alignment is broken -_- Perhaps you must choose the right
    %CS frames to align to. Add in the trial number 
   %
   %TODO PLOT AVERAGE WRITES AIR
    %USE WITH V6 ANALYSIS
    %{
    INPUTS:
    trials: The entire dataset segmented by trials, not including ITIS.
    pre_ITIS: The ITIs before each trial start.
    post_ITIs: The ITIs after each trial start.
    baseline_frames: The number of frames to use in calculating background.
    seconds_per_trial: ... seconds per trial
    alignment_frame: cell of frames where the alignment event is delivered.
        ex. If you use water_frame, that will align to water US events.
    CS_frame: The cell of frames containing CS onset info for each trial.
        It's used to calculate CS_target_frame.
    trial_by_trial: 1 to calculate baselines trial by trial. 0 to use an 
        estimated global baseline
    average_baseline_signal: The estimated global baseline. If the previous
        variable is 0, this input doesn't matter.
    
    OUTPUTS:
    uniform_trials: All trials, aligned to alignment_frame
    US_trials: All trials with a valid alignment frame. For water_frame, 
        this is trials that have water delivered. For sound_frame, this is
        every trial with sound delivered.
    no_US_trials: All trials without an alignment frame.
    US_target_frame, CS_target_frame
    %}
    
    %TODO write in an extend option instead of this greater than 10
    %business
%    no_US_trials = [];
%    trials_to_cut = [];
%    if strcmp(alignment, 'escape')
%        trials_to_cut = find(which_trials == 0);
%        which_trials = ones(1, sum(which_trials));
%    end
    
    %if ~isempty(trials_to_cut)
    %    trials(trials_to_cut)=[];
    %    post_ITIs(trials_to_cut) = [];
    %    pre_ITIs(trials_to_cut) = [];
    %    US_frame(trials_to_cut) = [];
    %    CS_frame(trials_to_cut) = [];
    %end    
    
    %You need to not consider the trials you are not looking at
    %US_frame{~which_trials} = 0;
    %for x = 1:size(which_trials, 2)
    %    if ~which_trials(x)
    %        US_frame{x} = 0;
    %        if iscell(CS_frame)
    %            CS_frame{x} = 0;
    %        else
    %            CS_frame(x) = 0;
    %        end
    %    end
    %end
    uniform_trials = [];
    sound_duration = 1;
    frame_offsets = [];
    true_sec_offsets = [];
    %Find the smallest trial size
    smallest_size = 100000;
    sec_smallest_size = 99999;
    trial_num = size(trials,2); 
    for s = 1:trial_num
       smoll = size(cell2mat(trials(1, s)), 1);
       if and(smallest_size > smoll, smallest_size ~= 0)
           sec_smallest_size = smallest_size;
           smallest_size = smoll;
       end
    end
   
    
    extension = 0;
    if seconds_to_extend > 0
       
       extra_frames = (seconds_to_extend) * frame_rate;
       extension = round(extra_frames);
    end
    
    %Change this to change the number of frames to extend 
                    %trial length by. Must be less than minimum ITI length.
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate the average frame in each trial that US/CS is received
    non_zero_US_frames = [];
    %if iscell(US_frame)
        %US_frame = cell2mat(US_frame) .* which_trials;
    %else
    %    US_frame = US_frame .* which_trials;
    %end
    for x = US_frame
        if iscell(US_frame)
            if x{1} ~= 0
                non_zero_US_frames = horzcat(non_zero_US_frames, x(1));
            end
        else
            if x(1) ~= 0
                non_zero_US_frames = horzcat(non_zero_US_frames, x(1));
            end
        end
    end
    
    if iscell(non_zero_US_frames)
        US_target_frame = round(sum(cell2mat(non_zero_US_frames))/length(non_zero_US_frames));
    else
        US_target_frame = round(sum(non_zero_US_frames)/length(non_zero_US_frames));
    end

    true_seconds_per_trial = seconds_per_trial + seconds_to_extend;%/smallest_size * (smallest_size + extension);
    
    %Set up the target alignment frame
    average_baselines = zeros(1, trial_num);
    average_signal = zeros(smallest_size, 1);

    alignment_frame = US_frame;
    target_frame = US_target_frame;
    
    %Set up matrices to hold trial information
    if iscell(US_frame)
        n_no_US_trials = sum(cell2mat(US_frame) == 0);
        US_frame = cell2mat(US_frame);
        %if ~strcmp(alignment, 'baseline run')
       %     US_frame = US_frame(find(which_trials == 1));
       % end
        
    else
        n_no_US_trials = sum(US_frame == 0);
        %if ~strcmp(alignment, 'baseline run')
        %    US_frame = US_frame(find(which_trials == 1));
        %end
    end
    n_US_trials = trial_num - n_no_US_trials;
    US_trials = zeros(n_US_trials, smallest_size + extension);
    uniform_trials = zeros(size(trials, 2), smallest_size + extension);
    %no_US_trials = zeros(n_no_US_trials, smallest_size + extension);

    %Make sure each trial is aligned correctly to bpod data
    reward_count = 1;
    fail_count = 1;
    uniform_ind = 1;
    for t = 1:trial_num
        single_trial = cell2mat(trials(t));
        pre = vertcat(zeros(120, 2), cell2mat(pre_ITIs(t)));
        post = vertcat(cell2mat(post_ITIs(t)), zeros(200, 2));
        whole_trial = vertcat(single_trial, post);
       
        if length(single_trial(:, 2)) > smallest_size %If too big 
        %Trim trials to the smallest size and calculate dF/F
            trial_signal = single_trial(1:smallest_size, :); %- average_baselines(t)' + auto_flo)/(average_baselines(t)' - auto_flo);
            pre_signal = pre; %- average_baselines(t)' + auto_flo)/(average_baselines(t)' - auto_flo);
            spill = single_trial((smallest_size + 1):end, :);
            post_signal = vertcat(spill, post); %- average_baselines(t)' + auto_flo)/(average_baselines(t)' - auto_flo);
            average_signal = average_signal + trial_signal(:,2);
        else  %If the correct size
        %calculate dF/F
            trial_signal = single_trial(:, :);% - average_baselines(t)' + auto_flo)/(average_baselines(t)' - auto_flo);
            pre_signal = pre;% - average_baselines(t)' + auto_flo)/(average_baselines(t)' - auto_flo);
            post_signal = post;% - average_baselines(t)' + auto_flo)/(average_baselines(t)' - auto_flo);
            average_signal = average_signal + trial_signal(:,2);
        end
        
        %Align each trial's onset of CS or US
        if iscell(US_frame)
            frame_U = US_frame{t};
        else
            frame_U = US_frame(t);
        end
        align_frame = alignment_frame(t);
        if iscell(align_frame)
            align_frame = cell2mat(align_frame);
        end
        %If trials are not aligned
        if and(and(align_frame ~= target_frame, align_frame ~= 0), 1)
            offset = align_frame - target_frame;
            lick_time = whole_trial(align_frame, 1);
            correct_time = whole_trial(target_frame, 1);
            true_sec_offsets = horzcat(true_sec_offsets, lick_time - correct_time);
            frame_offsets = horzcat(frame_offsets, offset);
            %If trial is too shifted right, add post ITI frames.
            if offset > 0 
                overshot = 0;
                if abs(offset) > size(trial_signal, 1) %If offset is too big do nothing
                    disp(horzcat('Offset is too big: ', num2str(offset)));
                    overshot = offset - size(trial_signal, 1);
                    if size(post_signal, 1) < (extension + overshot + size(trial_signal, 1))
                        continue
                    end
                end
                if seconds_to_extend > 0
                    if overshot ~= 0
                        halp = post_signal(overshot + 1:(extension + overshot + size(trial_signal, 1)), 2); removed size(trialsignal)
                        %halp = post_signal(overshot + 1:(extension + overshot), 2);
                        uniform_trials(uniform_ind,:) = halp;
                    else
                        if size(post_signal, 1) < (extension + offset)
                            error('choose a shorter seconds_to_extend');
                        end
                        uniform_trials(uniform_ind,:) = vertcat(trial_signal(offset + 1:end, 2), post_signal(1:offset, 2), post_signal((offset + 1):(extension + offset), 2));
                    end
                else
                    if offset >= smallest_size
                        diff = offset - smallest_size;
                        uniform_trials(uniform_ind,:) = post_signal((1 + diff):(smallest_size + diff), 2);
                    else
                        uniform_trials(uniform_ind,:) = vertcat(trial_signal(offset + 1:end, 2), post_signal(1:offset, 2));
                    end
                end
                uniform_ind = uniform_ind + 1;
            %If too shifted left, add pre ITI frames.
            else 
                if seconds_to_extend > 0
                    if abs(offset) > 20 %If offset is too big do nothing
                        disp(horzcat('Offset is too big: ', num2str(offset)));
                        %if seconds_to_extend > 0
                        %    uniform_trials(t,:) = vertcat(trial_signal, post_signal(1:extension));
                        %else
                        %    uniform_trials(t,:) = trial_signal;
                        %end
                        %continue;
                    end
                    a = pre_signal((end - abs(offset) + 1):end, 2); 
                    b = trial_signal(1:(end + offset), 2);          
                    c = trial_signal((end + offset + 1):end, 2); %Don't know why you'd need excess    
                    d = post_signal(1:(extension - abs(offset)), 2);
                    uniform_trials(uniform_ind,:) = vertcat(a, b, c, d);
                else
                    %%%%%%%%TODO FIX INDEXING ERROR
                    preprepre = transpose(pre_signal((end + offset + 1):end, 2));
                    postpostpost = transpose(trial_signal(1:(end + offset), 2));
                    
                    uniform_trials(uniform_ind,:) = horzcat(preprepre, postpostpost);
                end
                uniform_ind = uniform_ind + 1;
            end  
        %If the trial is already aligned
        else   %TODO if you want, keep track of those with no licks
            if seconds_to_extend > 0
                uniform_trials(uniform_ind,:) = vertcat(trial_signal(:,2), post_signal(1:extension, 2));
            else
                uniform_trials(uniform_ind,:) = trial_signal(:,2);
            end
            frame_offsets = horzcat(frame_offsets, 0);
            %lick_time = single_trial(align_frame, 1);
            %correct_time = single_trial(target_frame, 1);
            true_sec_offsets = horzcat(true_sec_offsets, 0);
            uniform_ind = uniform_ind + 1;
        end 
        %Split the trials by whether or not there was a US
        %if frame_U == 0  
        %    no_US_trials(fail_count,:) = uniform_trials(t,:);
        %    fail_count = fail_count + 1;
        %else
        %    US_trials(reward_count,:) = uniform_trials(t,:);
        %    reward_count = reward_count + 1;
        %end
    end
    if uniform_ind < 1 + size(uniform_trials(1))
        uniform_trials((uniform_ind:end), :) = [];
    end
    sec_offsets = frame_offsets ./ frame_rate;
end

