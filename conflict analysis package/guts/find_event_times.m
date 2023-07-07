function [event_times, event_end_times, event_trial_ind, no_event_trial_ind, event_labels] = find_event_times(SessionData, alt_start_times, alignment, trials_to_check, reward_amount, punish_amount, leave_out_no_licks, stim_time, sound_on, outcomes, left_lick_port, right_lick_port, lick_window)
    no_left_ind = [];
    left_ind = [];
    no_right_ind = [];
    right_ind = [];
    left_lick_event_times = [];
    right_lick_event_times = [];
    left_deliveries = [];
    right_deliveries = [];
    punish_i = 0;
    reward = 7;
    water_delivered_times = [];
    water_earned_times = [];
    water_earned_ind = [];
    no_water_delivered_ind = [];
    water_delivered_ind = [];
    water_assisted_times = [];
    water_assisted_ind = [];
    sound_delivered_times = [];
    sound_delivered_ind = [];
    air_sounds = [];
    air_sound_ind = [];
    water_sounds = [];
    water_sound_ind = [];
    avoid_times = [];
    avoid_ind = [];
    not_avoided = [];
    total_not_avoided = [];
    run_times = [];
    escape_times = [];
    escape_ind = [];
    laser_times = [];
    laser_ind = [];
    air_puff_times = [];
    air_puff_ind = [];
    no_air_delivered = [];
    air_puff_end_times = [];
    air_lengths = [];
    air_fail_ind = [];
    %air_trial_ind = 0;
    event_end_times = [];
    RawEvent = SessionData.RawEvents.Trial;
    first_volume = find_first_volume(SessionData);
    start_times = alt_start_times;
    no_event_trial_ind = [];
    event_trial_ind = [];
    licky = [];
    outcomes = fix_outcomes(SessionData);
    trial_types = [];
    if isfield(SessionData, 'TrialTypes')
        if ~isempty(SessionData.TrialTypes)
            trial_types = SessionData.TrialTypes;
        else
            trial_types = 1:SessionData.nTrials;
            trial_types = mod(trial_types, 2);
            outcomes = ~trial_types;
        end
    else
        spout = SessionData.spout;
        if spout == 3
            trial_types = 1:SessionData.nTrials;
            trial_types = mod(trial_types, 2);
            outcomes = ~trial_types;
        elseif spout == 2
            trial_types = ones(1, SessionData.nTrials);
            outcomes = ~trial_types;
        elseif spout == 1
            trial_types = zeros(1, SessionData.nTrials);
            outcomes = ~trial_types;
        end
    end
    %right_stim_time = SessionData.RawEvents.Trial{1, 2}.States.DeliverRight(2);
%    left_stim_time = SessionData.RawEvents.Trial{1, 2}.States.DeliverLeft(2);
    %It's not choosing the right kind of trial
    %3 receives reward with assist. 1 receive reward without assist. 0.
    %punishment. -1, avoided air puff. 2 for ignored reward. 4 escape 
    if ~strcmp(alignment, 'events')
        %trial_types = SessionData.TrialTypes;
        for q = trials_to_check
            trial_start = start_times(q);  %%%% CHANGE THIS to new start times
            if or(strcmp(alignment, 'both'), or(strcmp(alignment, 'left'), strcmp(alignment, 'right')))
                if isfield(RawEvent{q}.Events,left_lick_port)
                    left_lick_event_times = horzcat(left_lick_event_times, eval(horzcat('RawEvent{q}.Events.', left_lick_port)) + trial_start);
                    %licky = left_lick_event_times;
                end
                if isfield(RawEvent{q}.Events,right_lick_port)
                    right_lick_event_times = horzcat(right_lick_event_times, eval(horzcat('RawEvent{q}.Events.', right_lick_port)) + trial_start);
                    %licky = left_lick_event_times;
                end
            else
                if isfield(RawEvent{q}.Events,left_lick_port)
                    lick_event_times = horzcat(left_lick_event_times, eval(horzcat('RawEvent{q}.Events.', left_lick_port)) + trial_start);
                end
            end
            if strcmp(alignment, 'left')
                reward_state = 'DeliverLeft';
            elseif strcmp(alignment, 'right')
                reward_state = 'DeliverRight';
            elseif strcmp(alignment, 'both')
            else
                reward_state = 'DeliverReward';
            end
            %Record left delivery times
            if trial_types(q) == 0
                if or(strcmp(alignment, 'left'), strcmp(alignment, 'right'))
                    reward_state = 'DeliverLeft';
                elseif strcmp(alignment, 'both')
                    reward_state = 'DeliverLeft';
                end
                if isfield(SessionData.RawEvents.Trial{1, q}.States, reward_state)
                    if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', reward_state, '(1)')))
                        reward_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', reward_state, '(1)'));
                        if q ~= 1
                            if isfield(SessionData, 'TrialSettings')
                                if isfield(SessionData.TrialSettings(q - 1).GUI, 'RewardAmount')
                                    reward = SessionData.TrialSettings(q - 1).GUI.RewardAmount;
                                elseif isfield(SessionData.TrialSettings(q - 1).GUI, 'LeftAmount')
                                    reward = SessionData.TrialSettings(q - 1).GUI.LeftAmount;
                                else
                                    reward = first_volume;
                                end
                            else
                                if strcmp(reward_state, 'DeliverLeft')
                                    reward = SessionData.LeftAmount;
                                else
                                    reward = SessionData.RightAmount;
                                end
                            end
                        else
                            reward = first_volume;
                        end
                        if and(~isnan(reward_time), or(reward == reward_amount, reward_amount == 1000))
                            if leave_out_no_licks
                                %Find trials with no licking after the US
                                licky = left_lick_event_times - trial_start;
                                left_stim_time = SessionData.RawEvents.Trial{1, q}.States.DeliverLeft(1);
                                if isfield(SessionData.RawEvents.Trial{1, q}.States, 'LickCheck')
                                    if ~isnan(SessionData.RawEvents.Trial{1, q}.States.LickCheck(2))
                                        left_stim_time = SessionData.RawEvents.Trial{1, q}.States.LickCheck(2);
                                        lick_window = 5; 
                                    end
                                end
                                if ~isempty(licky)
                                    licky = licky(find(licky >= (left_stim_time))); %9
                                end
                                if ~isempty(licky)
                                   licky = licky(find(licky < left_stim_time + lick_window)); %10 changed from
                                end                
                                if ~isempty(licky)
                                    %Exclude trials with no licking
                                    disp(q);
                                    first_lick = licky(1);
                                    reward_time = first_lick;
                                    if outcomes(q) == 3
                                        water_assisted_times = horzcat(water_assisted_times, reward_time + trial_start);
                                        water_assisted_ind = horzcat(water_assisted_ind, q);
                                        water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                        water_delivered_ind = horzcat(water_delivered_ind, q);
                                        left_deliveries = horzcat(left_deliveries, reward_time + trial_start);
                                        left_ind = horzcat(left_ind, q);
                                    else
                                        water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                        left_deliveries = horzcat(left_deliveries, reward_time + trial_start);
                                        left_ind = horzcat(left_ind, q);
                                        water_delivered_ind = horzcat(water_delivered_ind, q);
                                        water_earned_times = horzcat(water_earned_times, reward_time + trial_start);
                                        water_earned_ind = horzcat(water_earned_ind, q);
                                    end
                                else
                                    disp(horzcat('DELETED FOR NO LICKING: ', num2str(q)));
                                end
                            else
                                if outcomes(q) == 3
                                    water_assisted_times = horzcat(water_assisted_times, reward_time + trial_start);
                                    water_assisted_ind = horzcat(water_assisted_ind, q);
                                    water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                    water_delivered_ind = horzcat(water_delivered_ind, q);
                                    left_deliveries = horzcat(left_deliveries, reward_time + trial_start);
                                    left_ind = horzcat(left_ind, q);
                                else
                                    water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                    water_delivered_ind = horzcat(water_delivered_ind, q);
                                    water_earned_times = horzcat(water_earned_times, reward_time + trial_start);
                                    water_earned_ind = horzcat(water_earned_ind, q);
                                    left_deliveries = horzcat(left_deliveries, reward_time + trial_start);
                                    left_ind = horzcat(left_ind, q);
                                end
                            end
                        end
                    else
                        no_water_delivered_ind = horzcat(no_water_delivered_ind, q);
                        %left_deliveries = horzcat(left_deliveries, reward_time + trial_start);
                        no_left_ind = horzcat(no_left_ind, q);
                    end
                else
                    no_water_delivered_ind = horzcat(no_water_delivered_ind, q);
                end

            else
                %If right or air
                reward_state = 'DeliverRight';
                if or(strcmp(alignment, 'left'), strcmp(alignment, 'right'))
                    reward_state = 'DeliverRight';
                elseif strcmp(alignment, 'both')
                    reward_state = 'DeliverRight';
                end

                if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', reward_state, '(1)')))
                    reward_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', reward_state, '(1)'));
                    if q ~= 1
                        if isfield(SessionData, 'TrialSettings')
                            if isfield(SessionData.TrialSettings(q - 1).GUI, 'RewardAmount')
                                reward = SessionData.TrialSettings(q - 1).GUI.RewardAmount;
                            elseif isfield(SessionData.TrialSettings(q - 1).GUI, 'RightAmount')
                                reward = SessionData.TrialSettings(q - 1).GUI.RightAmount;
                            else
                                reward = first_volume;
                            end
                        else
                            reward = 1;
                        end
                    else
                        reward = first_volume;
                    end
                    if and(~isnan(reward_time), or(reward == reward_amount, reward_amount == 1000))
                        if leave_out_no_licks
                            %Find trials with no licking after the US
                            licky = right_lick_event_times - trial_start;
                            right_stim_time = SessionData.RawEvents.Trial{1, q}.States.DeliverRight(1);
                            if isfield(SessionData.RawEvents.Trial{1, q}.States, 'LickCheck')
                                if ~isnan(SessionData.RawEvents.Trial{1, q}.States.LickCheck(2))
                                    right_stim_time = SessionData.RawEvents.Trial{1, q}.States.LickCheck(2);
                                    lick_window = 5; 
                                end
                            end
                            if ~isempty(licky)
                                licky = licky(find(licky >= (right_stim_time))); %9
                            end
                            if ~isempty(licky)
                               licky = licky(find(licky < right_stim_time + lick_window)); %10 changed from
                            end

                            if ~isempty(licky)
                                %Exclude trials with no licking
                                first_lick = licky(1);
                                reward_time = first_lick;
                                disp(q);
                                if outcomes(q) == 3
                                    water_assisted_times = horzcat(water_assisted_times, reward_time + trial_start);
                                    water_assisted_ind = horzcat(water_assisted_ind, q);
                                    water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                    water_delivered_ind = horzcat(water_delivered_ind, q);
                                    right_deliveries = horzcat(right_deliveries, reward_time + trial_start);
                                    right_ind = horzcat(right_ind, q);
                                else
                                    water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                    water_delivered_ind = horzcat(water_delivered_ind, q);
                                    water_earned_times = horzcat(water_earned_times, reward_time + trial_start);
                                    water_earned_ind = horzcat(water_earned_ind, q);
                                    right_deliveries = horzcat(right_deliveries, reward_time + trial_start);
                                    right_ind = horzcat(right_ind, q);
                                end
                            else
                                disp(horzcat('DELETED FOR NO LICKING: ', num2str(q)));
                            end
                        else
                            if outcomes(q) == 3
                                water_assisted_times = horzcat(water_assisted_times, reward_time + trial_start);
                                water_assisted_ind = horzcat(water_assisted_ind, q);
                                water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                water_delivered_ind = horzcat(water_delivered_ind, q);
                                right_deliveries = horzcat(right_deliveries, reward_time + trial_start);
                                right_ind = horzcat(right_ind, q);
                            else
                                water_delivered_times = horzcat(water_delivered_times, reward_time + trial_start);
                                water_delivered_ind = horzcat(water_delivered_ind, q);
                                water_earned_times = horzcat(water_earned_times, reward_time + trial_start);
                                water_earned_ind = horzcat(water_earned_ind, q);
                                right_deliveries = horzcat(right_deliveries, reward_time + trial_start);
                                right_ind = horzcat(right_ind, q);
                            end
                        end
                    end
                else
                    no_water_delivered_ind = horzcat(no_water_delivered_ind, q);
                    no_right_ind = horzcat(no_right_ind, q);
                end
            
            end
            %Record sound delivered times
            if isfield(SessionData.RawEvents.Trial{1, q}.States, 'StimulusDeliver')
                if ~isnan(SessionData.RawEvents.Trial{1, q}.States.StimulusDeliver(1))
                    stimulus_time = SessionData.RawEvents.Trial{1, q}.States.StimulusDeliver(1);
                    sound_delivered_times = horzcat(sound_delivered_times, stimulus_time + trial_start);% - session_start);
                    sound_delivered_ind = horzcat(sound_delivered_ind, q);% - session_start);
                    if trial_types(q) == 0 %water
                        water_sounds = [];
                        water_sound_ind = [];
                    elseif trial_types(q) == 1
                        air_sounds = [];
                        air_sound_ind = [];
                    end
                else
                    sound_delivered_times = horzcat(sound_delivered_times, 0); %Record trials with no sound as the start of the trial
                end
            end

            %Avoid time
            if isfield(SessionData.RawEvents.Trial{1, q}.States, 'DelayPre')
                if ~isnan(SessionData.RawEvents.Trial{1, q}.States.DelayPre(1))
                    avoid_time = SessionData.RawEvents.Trial{1, q}.States.DelayPre(1);
                    avoid_times = horzcat(avoid_times, avoid_time + trial_start);% - session_start);
                    avoid_ind = horzcat(avoid_ind, q);
                else
                   not_avoided = horzcat(not_avoided, q);
                   total_not_avoided = horzcat(total_not_avoided, q);
                end
            end

            %Running initiation times
            if isfield(SessionData.RawEvents.Trial{1, q}.Events, 'Wire1High')
                if ~isempty(SessionData.RawEvents.Trial{1, q}.Events.Wire1High) %Might need to change to isfield
                    run_time = SessionData.RawEvents.Trial{1, q}.Events.Wire1High;

                    run_times = horzcat(run_times, run_time + trial_start);% - session_start);

                end
            end

            %Escape time
            if isfield(SessionData.RawEvents.Trial{1, q}.States, 'Escape')
                if ~isnan(SessionData.RawEvents.Trial{1, q}.States.Escape(1))
                    escape_time = SessionData.RawEvents.Trial{1, q}.States.Escape(1);
                    escape_times = horzcat(escape_times, escape_time + trial_start);% - session_start);
                    escape_ind = horzcat(escape_ind, q);
                end
            end

            %Record laser times
            if isfield(SessionData.RawEvents.Trial{1, q}.States, 'Laser')
                if ~isnan(SessionData.RawEvents.Trial{1, q}.States.Laser(1))
                    stimulus_time = SessionData.RawEvents.Trial{1, q}.States.Laser(1);
                    laser_times = horzcat(laser_times, stimulus_time + trial_start);% - session_start);
                    laser_ind = horzcat(laser_ind, q);
                end
            end

            %Record Air puff times
            if trial_types(q) == 1
                if isfield(SessionData.RawEvents.Trial{1, q}.States, 'DeliverPunishment')
                    if ~isnan(SessionData.RawEvents.Trial{1, q}.States.DeliverPunishment(1))
                        %Keeps track of the air puff length
                        if q ~= 1
                            if isfield(SessionData, 'TrialSettings')
                                punish = SessionData.TrialSettings(q - 1).GUI.RewardAmount;
                            else
                                punish = 1;
                            end
                        else
                            punish = first_volume;
                        end
                        punishment_time = SessionData.RawEvents.Trial{1, q}.States.DeliverPunishment(1);
                        punishment_end = SessionData.RawEvents.Trial{1, q}.States.DeliverPunishment(2);
                        if and(~isnan(punishment_time), or(punish == punish_amount, punish_amount == 1000))
                            %record the air puff if length is right
                            air_puff_times = horzcat(air_puff_times, punishment_time + trial_start);
                            air_puff_ind = horzcat(air_puff_ind, q);
                            air_puff_end_times = horzcat(air_puff_end_times, punishment_end + trial_start);
                            air_lengths = horzcat(air_lengths, punishment_end - punishment_time);
                            air_fail_ind = horzcat(air_fail_ind, q);
                        end
                    
                    else
                        no_air_delivered = horzcat(no_air_delivered, q);
                    end
                else
                    no_air_delivered = horzcat(no_air_delivered, q);
                end
            end
            if trial_types(q)
                punish_i = punish_i + 1;
                %air_trial_ind = horzcat(air_trial_ind, q);
            end
            left_lick_event_times = [];
            right_lick_event_times = [];
        end
    else
        trial_types = ones(1, length(trials_to_check)) .* 4;
    end

    %Find running initiation times. This is any running bought at least
    %1 second after the last running bought and continues for 3 or .
    run = -1;
    initiations = [];
    bought_length = 3;
    bought_lengths = [];
    for p = 1:size(run_times,2)
        next_run = run_times(1,p);
        if next_run - run > 1
            if bought_length < 3
                %If the bought length was too short, delete that run
                if ~isempty(initiations)
                    initiations = initiations(1:end-1);
                end
            else
                if or(~isempty(initiations), p == size(run_times,2))
                    bought_lengths = horzcat(bought_lengths, bought_length);
                end
            end
            %Add the next running initiation
            initiations = horzcat(initiations, next_run);
            %reset bought length to 1
            bought_length = 1;
        elseif p == size(run_times,2)
            bought_lengths = horzcat(bought_lengths, bought_length);
        elseif next_run - run < .1
            bought_length = bought_length + 1;
            %Next make it so at the start of the next iteration, if
            %bought length is less than 2, dlete previous bought.
        elseif next_run - run < 1
            disp('I dont know what to do here');

        end
        run = next_run;
    end

    %Cycle through start times. Segment from start_times(q) to q + 1.
    %include everything < 4 seconds as baseline.
    first_start = start_times(1);
    baseline_initiations = [];
    for q = 2:size(start_times, 2)
        next_start = start_times(q);
        %Find all running initiations in the current trial
        trial_segment = initiations(find(and(initiations > first_start, initiations < next_start)));
        %Put into time scale of 10 second trial
        trial_segment = trial_segment - first_start;
        %Find all places where running initiates during baseline
        trial_segment = trial_segment(find(trial_segment < 4)); %trial_segment(find(or(trial_segment < 4, trial_segment > 10)));
        %Convert back to timescale of imaging data
        trial_segment = trial_segment + first_start;
        %Keep track of the initiations
        baseline_initiations = horzcat(baseline_initiations, trial_segment);
        first_start = next_start;
    end
    if sound_on
        event_labels = {'Signal', 'CS', alignment};
    else
        event_labels = {'Signal', alignment};
    end
    
    if strcmp(alignment, 'events')
        [LeftPokes, MidPokes, RightPokes] = ExtractPokeTimes(SessionData);
        which_event = input('Left (l), middle (m), or right (r)? ', 's');
        if strcmp(which_event, 'l')
            event_times = {LeftPokes};
            event_labels = {'Signal', foods{1}};
        elseif strcmp(which_event, 'm')
            event_times = {MidPokes};
            event_labels = {'Signal', foods{2}};
        elseif strcmp(which_event, 'r')
            event_times = {RightPokes};
            event_labels = {'Signal', foods{3}};
        end
        event_trial_ind = 1:size(event_times{1}, 2);
    elseif strcmp(alignment, 'water all')
        event_times = {water_delivered_times};
        no_event_trial_ind = no_water_delivered_ind;
        event_trial_ind = water_delivered_ind;
    elseif strcmp(alignment, 'assist')
        event_times = {water_assisted_times};
        no_event_trial_ind = no_water_delivered_ind;
        event_trial_ind = water_assisted_ind;
    elseif strcmp(alignment, 'water')
        event_times = {water_earned_times};
        no_event_trial_ind = no_water_delivered_ind;
        event_trial_ind = water_earned_ind;
    elseif strcmp(alignment, 'air')
        event_times = {air_puff_times};
        event_end_times = air_puff_end_times;
        no_event_trial_ind = no_air_delivered;
        event_trial_ind = air_puff_ind;
    elseif strcmp(alignment, 'baseline run')
        event_times = {baseline_initiations};
    elseif strcmp(alignment, 'escape')
        event_times = {escape_times};
        event_trial_ind = escape_ind;
    elseif strcmp(alignment, 'avoid')
        event_times = {avoid_times};
        event_trial_ind = avoid_ind;
    elseif strcmp(alignment, 'laser')
        event_times = {laser_times};
        event_trial_ind = laser_ind;
    elseif strcmp(alignment, 'air sound')
        event_times = air_sounds;
        event_trial_ind = air_sound_ind;
    elseif strcmp(alignment, 'left')
        event_times = {};
        event_times{1} = left_deliveries;
        event_times{2} = right_deliveries;
        event_trial_ind = left_ind;
        %no_event_trial_ind = right_ind;
        no_event_trial_ind = setdiff(trials_to_check, left_ind);
    elseif strcmp(alignment, 'right')
        event_times = {};
        event_times{1} = right_deliveries;
        event_times{2} = left_deliveries;
        event_trial_ind = right_ind;
        %no_event_trial_ind = left_ind;
        no_event_trial_ind = setdiff(trials_to_check, right_ind);
    elseif strcmp(alignment, 'both')
        event_times = {};
        total_deliveries = horzcat(left_deliveries, right_deliveries);
        event_times{1} = sort(total_deliveries);
        event_times{2} = [];
        total_ind = horzcat(left_ind, right_ind);
        event_trial_ind = sort(total_ind);
        %no_event_trial_ind = left_ind;
        no_event_trial_ind = setdiff(trials_to_check, total_ind);
    end
    
    if sound_on
        event_times{2} = sound_delivered_times;
    end
    if trials_to_check(1) > 1
        %event_trial_ind = event_trial_ind - trials_to_check(1) + 1;
        %no_event_trial_ind = no_event_trial_ind - trials_to_check(1) + 1;
    end
end