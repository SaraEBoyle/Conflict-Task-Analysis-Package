%% Instructions:
%% Update: Fixed right spout lps. Added lps per trial, licks per trial. 
%saved in 'Behavior Figures' folder. Z score peaks are in 'peak graphs' 
%folder, 'first_peaks' and 'second_peaks' variables. I will make this more 
%intuitive later. 

function  Photometry_Analysis_V9()
    %% Clear Workspace
    %clear, clc
    for signal = 1:2
    close all
    to_analyze = {'left', 'right'}; %choose sound, water, light, shock
    nose_poke_port = 'Port3In';
    pre_time = 10;
    post_time = 10;
    seconds_per_trial = pre_time + post_time;
    seconds_to_extend = 0;
    trials_to_check = 2:50; %List of trials (Don't worry for event alignment)
    fiber = 1; %Choose which fiber to look at
    pre_ITI_length = 15;
    trigger1 = 1; %Change to 1 if your signal alternates through 2 channels
    frame_rate = 20.0;
    plot_single_trials = 0;
   
    %Set these so the program can identify your files
    signal_key_word = 'SIGNAL'; %Use a word that's in your signal files
    analog_key_word = 'ANALOG'; %Use a word in your analog files
    bpod_key_word = 'CSand'; %Use a word in the bpod files
      
    %% Set up Variables
    %What type of bleach correction do you want? Choose only 1
    moving_avg = 1; %Does a moving average to subtract F.
    correct_for_bleaching = 0; %Subtracts biexponential curve to correc
    
    %Variables to select what you want to do
    stim_duration = .2;
    subtract_autofluorescence = 1; % if 0, don't subtract. If 1, choose the values to subtract carefully
    plot_heatmap_and_average = 1; %Haven't worked this out when water_align variable is on 
 
    %Never change below here
    if signal == 1
        if volume_to_check == 1000
            fiber_name = 'signal';
        else
            fiber_name = horzcat(num2str(volume_to_check), ' uL signal');
        end
    else
        if volume_to_check == 1000
            fiber_name = 'control';
        else
            fiber_name = horzcat(num2str(volume_to_check), ' uL control');
        end
    end
    
    if volume_to_check == 1000
        vol_name = '';
    else
        vol_name = horzcat(num2str(volume_to_check), ' ');
    end
                
    %% Before the for loop   
    file = pwd;
    file = horzcat(file, '/');
    
    if trigger1
        frame_rate = frame_rate/2;
    end
    
    %% Load data
    [SessionData, data, start_times, alt_start_times] = load_data(file, bpod_key_word, signal_key_word, analog_key_word, 1, 1, 1);
    
    
    
    
    
    
    
    
    [data, auto_flo] = trim_data(data, fiber, subtract_autofluorescence);
    
    %De-interleave and subtract autofluorescence
    [gcamp, iso] = de_interleave(trigger1, 0, file, fiber_name, data, fiber, auto_flo, subtract_autofluorescence, signal);
    
    %correct bleaching
    [gcamp, iso] = correct_bleaching(correct_for_bleaching, gcamp, iso, moving_avg, frame_rate, trigger1, 0);
    
    %Is sound on?
    sound_on = 0;
    sound_only = 0;
    man_target_frame = seconds_per_trial + seconds_to_extend;
    man_target_frame = frame_rate * man_target_frame/2;
    man_target_frame = 100;
    for iteration = 1:size(to_analyze, 2)

        %% Don't change anything below here
        alignment = to_analyze{iteration};
        align_to_US = 1;
  
        if or(strcmp(alignment,'left'), strcmp(alignment,'right'))
            %To look at water trials
            water_on = 1; 
            air_on = 0;
        elseif strcmp(alignment,'both')
            water_on = 1; 
            air_on = 0;
        end  
        %Variables for saving data
        save_graphs = 1; 
        save_data = 1;
        create_save_folder(fiber_name, alignment, file, right_spout, left_spout, vol_name);
        mode = horzcat(num2str(sound_on), num2str(water_on), num2str(air_on));
        %Record time of each stimulus
        %if ~isnan(SessionData.Outcomes(1))
            %outcomes = SessionData.Outcomes;
        %else
            %outcomes = ones(1, size(trial_types, 2));
        %end
        if isfield(SessionData, 'TrialTypes')
        if ~isempty(SessionData.TrialTypes)
            trial_types = SessionData.TrialTypes;
        else
            trial_types = 1:SessionData.nTrials;
            trial_types = mod(trial_types, 2);
            outcomes = trial_types;
        end
        elseif isfield(SessionData, 'spout')
            spout = SessionData.spout;
            if spout == 3
                trial_types = 1:SessionData.nTrials;
                trial_types = mod(trial_types, 2);
                outcomes = trial_types;
            elseif spout == 2
                trial_types = ones(1, SessionData.nTrials);
            elseif spout == 1
                trial_types = zeros(1, SessionData.nTrials);
            end
        end
    
        outcomes = ones(1, size(trial_types, 2)); %added
        
        [event_times, event_end_times, event_trial_ind, no_event_trial_ind, event_labels] = find_event_times(SessionData, alt_start_times, alignment, trials_to_check, volume_to_check, volume_to_check, leave_out_no_licks, stim_time, sound_on, outcomes, nose_poke_port, right_lick_port, lick_window);
         if isempty(event_trial_ind)
             disp(horzcat(alignment, ' ', num2str(volume_to_check),' uL trials do not exist'))
             cd(file);
             if strcmp(alignment, 'left')
                 rmdir(horzcat(file, fiber_name, '/left ', left_spout));
             else
                 rmdir(horzcat(file, fiber_name, '/right ', right_spout));
             end
             continue
         end
        
        %% De-interleave and subtract autofluorescence      
        gcamp_clone = gcamp;
         
        %% Dividing trials based on start times
        if iscell(event_times)
            size_to_check = size(event_times, 2);
        else
            size_to_check = size(event_times, 1);
        end
    
        if and(size_to_check > 1, ~strcmp(alignment, 'right'))
            sound_delivered_times = event_times{2};
            event_times = event_times{1};
        else
            sound_delivered_times = [];
            event_times = event_times{1};
        end
        %for all trials, not just US receiving
        trials = {};
        iso_trials = {};
        OG_trials = {};
        pre_ITIs = {};
        post_ITIs = {};
        US_frame = {}; 
        no_US_delivered = [];
        US_delivered = [];
        ind = 1;
        US_ind = 1;
        US_frames = [];
        CS_frames = [];
        checked_trial_types = [];
        other_side = [];
        %Just for US receiving
        event_trials = {};
        event_i = 1;
        event_pre_ITIs = {};
        event_post_ITIs = {};
        event_times_temp = event_times;
        event_frame = {};          %Record frame air puff happens                              
        event_delivered = [];
        
        %Just for no US trials
        no_event_trials = {};
        no_event_i = 1;
        no_event_pre_ITIs = {};
        no_event_post_ITIs = {};
        sound_delivered_times_temp = sound_delivered_times;
        %start_times = alt_start_times; 
        start_times = event_times - pre_time;
        end_times = event_times + post_time;
        evs = 1;
        %TODO INVESTIGATE
        for y = 1:length(start_times)
            %Segment the trials based on start and end times
            %pre_ITI_length = 5;
            post_ITI_length = 150;   
            if evs > length(event_times)
                continue
            end
            starty = start_times(y);
            endy = end_times(y);
            indices1 = find(gcamp(:,1) >= starty);
            if trigger1
                indices2 = find(iso(:,1) >= starty);
                trial2= iso(indices2,:);
            end
            indices3 = find(gcamp_clone(:,1) >= starty);
            trial1=gcamp(indices1,:);
            trial3 = gcamp_clone(indices3, :);
            
            indices1 = find(trial1(:,1) <= endy);
            if trigger1
                indices2 = find(trial2(:,1) <= endy);
                trial2=trial2(indices2,:);
            end
            indices3 = find(trial3(:,1) <= endy);
            trial1=trial1(indices1,:);
            trial3=trial3(indices3,:);
            pre_ITI=gcamp(find(gcamp(:,1) >= (start_times(y) - pre_ITI_length)),:); %Bigger than start - 15
            pre_ITI=pre_ITI(find(pre_ITI(:,1) <= (start_times(y))), :); %smaller than start of trial?
            post_ITI= gcamp(find(gcamp(:,1) >= end_times(y)), :); %Greater than end
            post_ITI_clone= post_ITI(find(post_ITI(:,1) <= (end_times(y) + post_ITI_length)),:);  %Less than end + ITI
            post_ITI= post_ITI(find(post_ITI(:,1) <= (end_times(y) + post_ITI_length)),:);  %Less than end + ITI
            
            US_ind = event_i;
            if ~isempty(event_times_temp)            %If not too many USs
                if y ~= 1
                    if  and(event_times_temp(1) < starty, event_times_temp(1) > end_times(y - 1) + seconds_to_extend)
                        disp(horzcat('The mouse licked too late into trial ', num2str(y - 1), '. Skipping that event'));
                        event_times_temp = event_times_temp(2:end);
                        to_take_out = find(event_trial_ind == y - 1);
                        event_trial_ind(to_take_out:end) = event_trial_ind(to_take_out:end) - 1;
                        event_trial_ind(to_take_out) = [];
                        ind = ind - 1;
                    end
                elseif event_times_temp(1) > end_times(1) + seconds_to_extend
                    disp(horzcat('The mouse licked too late into trial ', num2str(y - 1), '. Skipping that event'));
                end
                US_frames = find(trial1(:, 1) >= event_times_temp(1)); %Find US frames
                 if isempty(US_frames)
                    US_frames = find(post_ITI_clone(:, 1) >= event_times_temp(1)); %Find US frames
                    US_frames = US_frames + size(trial1, 1);
                    if ~isempty(US_frames)
                        if US_frames(1) > size(trial1,1) + seconds_to_extend*10 -20
                            disp(horzcat('The mouse licked too late into trial ', num2str(event_trial_ind(y)), '. Skipping that event'));
                            event_times_temp = event_times_temp(2:end);
                            to_take_out = find(event_trial_ind == event_trial_ind(y));
                            event_trial_ind(to_take_out:end) = event_trial_ind(to_take_out:end) - 1;
                            event_trial_ind(to_take_out) = [];
                            ind = ind - 1;
                            US_frames = find(trial1(:, 1) >= event_times_temp(1)); %Find US frames
                            if isempty(US_frames)
                                disp('This mouse licked too late into previous trial then too late into current')
                            end
                        end
                    end
                end
            end
                  
            if and(sound_on == 1, ~isempty(sound_delivered_times_temp))
                %Find CS_frames
                sound_time = sound_delivered_times_temp(1);
                if sound_time ~= 0
                    CS_frames = find(trial1(:, 1) >= sound_time);
                    sound_delivered_times_temp = sound_delivered_times_temp(2:end); %Delete CS
                else
                    sound_delivered_times_temp = sound_delivered_times_temp(2:end); %Delete CS
                end
            end

            if ~isempty(CS_frames)                               %If there's a CS
                sound_frame{ind} = CS_frames(1,1);
            else
                sound_frame{ind} = 0;
            end
            CS_frames = [];

            if or(~isempty(US_frames), sound_only)                               %If there's a US
                if sound_only
                    US_frame{ind} = sound_frame{ind};
                else
                    if US_frames(1, 1) ~= 1
                        US_frame{ind} = US_frames(1, 1);                 %Find the US onset frame
                    else
                        US_frame{ind} = 0;
                    end
                end                
                US_delivered = horzcat(US_delivered, event_trial_ind(y));         %record trial # with US
                trials{ind} = trial1(:, [1 2]);                  %Record trial data in total list
                if trigger1
                    iso_trials{ind} = trial2(:, [1 2]);
                end
                OG_trials{ind} = trial3(:, [1 2]);
                pre_ITIs{ind} = pre_ITI;                         %Record pre ITI data
                post_ITIs{ind} = post_ITI;                       %Record post ITI data
                            
                if ismember(event_trial_ind(y), US_delivered)
                    event_trials{event_i} = trial1(:, [1 2]);        %Record trial as air trial
                    event_pre_ITIs{event_i} = pre_ITI;               %Record air pre ITI
                    event_post_ITIs{event_i} = post_ITI;             %Record are post ITI
                    event_frame{event_i} = US_frame{ind};          %Record frame air puff happens                              
                    event_delivered = horzcat(event_delivered, event_i);
                    event_i = event_i + 1;
                    event_times_temp = event_times_temp(2:end);
                
                end
                US_frames = [];                                  %Clear US frames
                ind = ind + 1;
            else                                                 %If there's no US
                US_frame{ind} = 0;
                other_side = horzcat(other_side, event_trial_ind(y));
                no_US_delivered = horzcat(no_US_delivered, ind);
                
                if ismember(event_trial_ind(y), no_event_trial_ind)
                    no_event_trials{no_event_i} = trial1(:, [1 2]);        %Record trial as air trial
                    no_event_pre_ITIs{event_i} = pre_ITI;               %Record air pre ITI
                    no_event_post_ITIs{event_i} = post_ITI;             %Record are post ITI
                    no_event_i = no_event_i + 1;
                end
                trials{ind} = trial1(:, [1 2]);
                if trigger1
                    iso_trials{ind} = trial2(:, [1 2]);
                end
                OG_trials{ind} = trial3(:, [1 2]);
                checked_trial_types = horzcat(checked_trial_types, trial_types(event_trial_ind(y)));
                pre_ITIs{ind} = pre_ITI;
                post_ITIs{ind} = post_ITI;
                ind = ind + 1;
            end
        end
        
        stimulus_times = event_times;
        check_trials = trials;
        if sound_on
            stimulus_times_temp = {sound_delivered_times(event_trial_ind), event_times};
            
        else
            stimulus_times_temp = {zeros(1, size(event_times, 2)), event_times};
        end
        if sound_on
            if ~isempty(no_event_trial_ind)
                no_stimulus_times_temp = {sound_delivered_times(no_event_trial_ind), []};
            end
        else
            no_stimulus_times_temp = {zeros(1, size(event_times, 2)), event_times};
        end
        
        %% Plot single trials 
        if plot_single_trials
            if ~isempty(event_trial_ind)
                plot_single_trials_V2(1, event_trials, alignment, event_labels, stimulus_times_temp, stim_duration);
            end
        end
        if plot_single_trials
            if ~isempty(no_event_trial_ind)
                plot_single_trials_V2(0, no_event_trials, alignment, event_labels, no_stimulus_times_temp, stim_duration);
            end
        end
        
        %% Sort and align trials.
        if plot_heatmap_and_average
            align_to_this = alignment;
            if ~isempty(event_trial_ind)
                good_pre_ITIs = pre_ITIs;
                good_post_ITIs = post_ITIs;
                [US_trials, US_target_frame, CS_target_frame, end_CS_frame, true_seconds_per_trial, sec_offsets, true_sec_offsets] = sort_and_align_trials_V7(event_trials, good_pre_ITIs, good_post_ITIs, ...
                    seconds_per_trial, US_frame, sound_frame, align_to_this, align_to_US, frame_rate, seconds_to_extend, stim_duration, 1, pre_time*10);
                frame_offsets = sec_offsets * frame_rate;
                if volume_to_check == 1000
                    vol_name = '';
                else
                    vol_name = horzcat(num2str(volume_to_check), ' ');
                end
                trial_start_times = alt_start_times;
                trial_order = trials_to_check;
                save(horzcat(file, vol_name, alignment, ' alignment data.mat'), 'sec_offsets', 'frame_offsets', 'event_trial_ind', 'no_event_trial_ind', 'true_sec_offsets', 'event_times', 'trial_start_times', 'trial_order');
            end
            no_US_trials = [];
            if ~isempty(event_trial_ind)
                US_trials = US_trials * 100;
                baselines = [];
                for r = 1:size(US_trials, 1)
                    cur = US_trials(r, :);
                    baseline = cur(10:((pre_time - 1)*frame_rate - 5));
                    baselines(r, :) = baseline;
                end
                plot_average(mode, fiber_name, US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                ylabel('dF/F (%)', 'FontSize', 20);
                if strcmp(alignment, 'left')
                    title(horzcat(fiber_name, ': % dF/F in ', left_spout, ' trials'));
                elseif strcmp(alignment, 'right')
                    title(horzcat(fiber_name, ': ', right_spout, ' trials'));
                else
                    title(horzcat(fiber_name, ': % dF/F in US delivered ', alignment, ' trials'));
                end

                if save_graphs
                    if strcmp(alignment , 'left')
                        savefig('left_average_dFF.fig')
                    elseif strcmp(alignment , 'right')
                        savefig('right_average_dFF.fig')
                    else
                        savefig('US_average.fig');
                    end
                end
            else
                disp(horzcat('Found no ', alignment, ' trials. Exited.'))
                cd(file);
                continue
            end
            
            %% Plot trial by trial responses
            %if and(~isempty(no_event_trial_ind), ~strcmp(alignment, 'assist'))
                %baseline_no_US_trials = no_US_trials(:,20:(US_target_frame - 20));
            %else
                baseline_no_US_trials = [];
                las_baseline_no_US_trials = [];
                no_las_baseline_no_US_trials = [];
                z_no_US_trials = [];
            %end
            if ~isempty(event_trial_ind)            
                plot_heatmap(US_trials, fiber_name, true_seconds_per_trial, US_target_frame, CS_target_frame, alignment);
                title(horzcat(fiber_name, ': ', alignment, ' trials'));
            
                if strcmp(alignment, 'left')
                    savefig('left_delivery_heatmap.fig');
                elseif strcmp(alignment, 'right')
                    savefig('right_delivery_heatmap.fig');
                else
                    savefig('US_heatmap.fig');
                end              
            end
           
            increments = length(US_trials);
            if increments == 0
                increments = length(no_US_trials);
            end
            x_axis = ((1:increments)/increments)*true_seconds_per_trial;
            if volume_to_check == 1000
                if align_to_US
                    laser_time = x_axis(US_target_frame);
                else
                    laser_time = x_axis(CS_target_frame);
                end
                x_axis = x_axis - laser_time;
                if ~isempty(event_trial_ind)
                    if strcmp(alignment, 'left')
                        left_trials = US_trials;
                        save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                    elseif strcmp(alignment, 'right')
                        right_trials = US_trials;
                        save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                    else
                        save('dFF_data.mat', 'US_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                    end
                end

                %US
                if isempty(event_trial_ind)
                    baseline_US_trials = [];
                else
                    for x = 1:size(US_trials, 2)
                        if isnan(US_trials(1,x))
                            US_trials(1,x) = 0;
                        end
                    end
                    baseline_US_trials = US_trials(:,20:(US_target_frame - 20));
                end
                if and(~isempty(no_event_trial_ind), ~isempty(no_US_trials))
                    %baseline_trials = vertcat(baseline_US_trials, baseline_no_US_trials);
                    baseline_trials = baselines;
                    std_avg = mean(std(baseline_trials));
                    avg = mean(mean(baseline_trials));
                    diffs = (no_US_trials - avg);
                    z_no_US_trials = (diffs ./ std_avg);
                    diffs = (US_trials - avg);
                    z_US_trials = (diffs ./ std_avg); 
                    plot_average(mode, fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                    title(horzcat(fiber_name, ': ', alignment, ' trials'));
                    ylabel('z-score(dF/F)', 'FontSize', 20);
                    if strcmp(alignment, 'left')
                        savefig('left_z_score_average.fig');
                    elseif strcmp(alignment, 'right')
                        savefig('right_z_score_average.fig');
                    else
                        savefig('US_z_score_average.fig');
                    end
                    if and(~isempty(no_event_trial_ind), ~isempty(no_US_trials))
                        plot_average(mode, fiber_name, z_no_US_trials, true_seconds_per_trial, CS_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                        title(horzcat(fiber_name, ': ', alignment, ' trials'));
                        ylabel('z-score(dF/F)', 'FontSize', 20);
                        if strcmp(alignment, 'left')
                            z_right_trials = z_no_US_trials;
                            z_left_trials = z_US_trials;
                            save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');                                 
                        elseif strcmp(alignment, 'right')
                            z_right_trials = z_US_trials;
                            z_left_trials = z_no_US_trials;
                            save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                        else
                            save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis', 'US_target_frame');
                        end
                    end
                else
                    %baseline_trials = vertcat(baseline_US_trials);
                    baseline_trials = baselines;
                    std_avg = mean(std(baseline_trials));
                    avg = mean(mean(baseline_trials));
                    diffs = (US_trials - avg);
                    z_US_trials = (diffs ./ std_avg); 
                    %post_laser = round(mean(end_CS_frame)) + 1;
                    %laser_effect = mean(z_US_trials(:,post_laser:post_laser + 20));
%                    save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis');

                    plot_average(mode, fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                    title(horzcat(fiber_name, ': ', alignment, ' trials'));
                    ylabel('z-score(dF/F)', 'FontSize', 20);
                    if strcmp(alignment, 'left')
                        savefig('left_z_score_average.fig');
                        save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                        z_left_trials = z_US_trials;
                        save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                    elseif strcmp(alignment, 'right')
                        savefig('right_z_score_average.fig');
                        save('dFF_data.mat', 'right_trials', 'US_target_frame', 'x_axis', 'baselines');
                        z_right_trials = z_US_trials;
                        save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                    else
                        savefig('US_z_score_average.fig');
                    end
                    if and(~isempty(no_event_trial_ind), ~isempty(no_US_trials))
                        plot_average(mode, fiber_name, z_no_US_trials, true_seconds_per_trial, CS_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                        title(horzcat(fiber_name, ': ', alignment, ' trials'));
                        ylabel('z-score(dF/F)', 'FontSize', 20);
                        if strcmp(alignment, 'left') 
                            save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                            z_right_trials = z_no_US_trials;
                            z_left_trials = z_US_trials;
                            save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                        elseif strcmp(alignment, 'right')
                            save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                            z_right_trials = z_US_trials;
                            z_left_trials = z_no_US_trials;
                            save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                        else
                            save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis', 'US_target_frame');
                        end
                    end
                end    
            else
                if align_to_US
                    if ~isnan(US_target_frame)
                        laser_time = x_axis(US_target_frame);
                    else
                        laser_time = x_axis(72);
                        US_target_frame = 72;
                    end
                else
                    laser_time = x_axis(CS_target_frame);
                end
                if strcmp(alignment, 'right')
                    right_trials = US_trials;
                elseif strcmp(alignment, 'left')
                    left_trials = US_trials;
                end
                %x_axis = x_axis - laser_time;
                laser_time = x_axis(US_target_frame);
                x_axis = x_axis - laser_time;
                if strcmp(alignment, 'right')
                    save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                elseif strcmp(alignment, 'left')
                    save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                else
                    save('dFF_data.mat', 'US_trials', 'no_US_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines');
                end
                
                %US
                if or(strcmp(alignment, 'right'), strcmp(alignment, 'left'))
                    %baseline_no_US_trials = no_US_trials(:,1:(US_target_frame - 20));
                    baseline_US_trials = US_trials(:,1:(US_target_frame - 20));
                end
                %baseline_trials = vertcat(baseline_US_trials, baseline_no_US_trials);
                baseline_trials = baselines;
                std_avg = mean(std(baseline_trials));
                avg = mean(mean(baseline_trials));
                %diffs = (no_US_trials - avg);
                %z_no_US_trials = (diffs ./ std_avg);
                diffs = (US_trials - avg);
                z_US_trials = (diffs ./ std_avg); 
                if ~isempty(event_trial_ind)
                    plot_average(mode, fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                    if strcmp(alignment, 'left')
                        correct_title = 'left';
                    elseif strcmp(alignment, 'right')
                        correct_title = 'right';
                    else
                        correct_title = alignment;
                    end

                    title(horzcat(fiber_name, ': ', correct_title, ' trials'));
                    ylabel('z-score(dF/F)', 'FontSize', 20);
                    if save_graphs
                        if strcmp(alignment, 'left')
                                savefig('left_z_score_average.fig');
                                z_left_trials = z_US_trials;
                        elseif strcmp(alignment, 'right')
                            savefig('right_z_score_average.fig');
                            z_right_trials = z_US_trials;
                        else
                            savefig('US_z_score_average.fig');
                        end
                    end
                end

                if strcmp(alignment, 'left')
                    save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                elseif strcmp(alignment, 'right')
                    save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                else
                    save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis', 'US_target_frame');
                end
            end
        end
        % It takes 0.0001 s for arduino to read analog input. .00025 s for bpod
        % state changes. .0001 to do output action
        cd(file);
    end
    random_message(2);
    end
end

%function [] = create_save_folder(fiber_name, alignment, file)

%function [data, auto_flo] = trim_data(data, fiber, subtract_autofluorescence)

%function [event_times, event_end_times, event_trial_ind, no_event_trial_ind, event_labels] = find_event_times(SessionData, alt_start_times, alignment, trials_to_check, reward_amount, punish_amount, leave_out_no_licks, stim_time, sound_on, outcomes, left_lick_port, right_lick_port, lick_window)

%function [gcamp, iso] = de_interleave(trigger1, trigger2, file, fiber_name, data, fiber, auto_flo, subtract_autofluorescence, signal_v)

%function [gcamp, iso] = correct_bleaching(correct_for_bleaching, gcamp, iso, moving_avg, frame_rate, trigger1, smoot)
