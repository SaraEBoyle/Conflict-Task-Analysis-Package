function [] = full_photometry_conflict_analysis()
    %% KEY UPDATES- can now select multiple folders to analyze at once. You
    %% don't have to select the background image anymore, select folders instead.
    %% As always, make sure you use the key words during daq (or add them to the 
    %% correct files), set the events you want to plot, and make sure you set your 
    %% fibers and fiber names.

    clearvars -except BpodSystem
    
    bpod_key_word = 'Conflict';
    analog_key_word = 'Time';
    signal_key_word = 'photometry';
    background_key_word = 'back';
    %                   1        2        3        4        5        6
    %event_labels = {'water', 'light', 'sound', 'laser', 'shock', 'shock received', 
    %      7                  8                          9                        10
    %'shock avoided', 'sound in shock received', 'sound in shock avoided', 'platform entries', 
    %    11                        12                     13                         14
    %'platform exits', 'sound platform entries', 'first platform entries', 'first platform exits after shock', 
    %      15                     16
    %'ITI platform entries', 'first water'};
    %events_to_analyze = [14]; %1:13; everything % [3, 8, 9] all sound; [1, 2, 4, 5, 6, 7, 10, 11, 12, 13] everything but sound, [10, 11, 12, 13] platform entries
    
    %Look at all possible alignments
    %basic_stimuli = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
    basic_stimuli = [1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
    minimal_stimuli = [1, 3, 11, 12, 13, 15, 16];
    just_shock = [6, 7];

    %Look at sound alignments for full sound duration
    full_sound = [3, 8, 9];

    %event_groups = {basic_stimuli, full_sound};

    event_groups = {basic_stimuli, full_sound};
    %event_groups = {5, 6};
    %Time before stimuli to check
    pre_times = {5, 20};

    %Time after stimuli to check
    post_times = {5, 40};

    %which timepoints to consider baseline
    baseline_times = {5, 5}; 

    %Names of fibers you're checking. Order matters here
    fiber_names = {'CeL Left', 'CeL Right'};
    %fiber_names = {'CeL', 'TS'};
    %List fibers you want to check in the analysis
    %signal_fibers = [1, 2];
    signal_fibers = [1, 2];
    %Your control cha1nnel- if using isosbestic channel it will be the same
    %as your signal fiber
    control_fiber = [1, 2];
    
    %% DON'T CHANGE BELOW HERE
    og_file = pwd;
    folders_to_analyze = uigetdir2(pwd);
    if isempty(folders_to_analyze)
        folders_to_analyze = {og_file};
    end
    for fold = 1:size(folders_to_analyze, 2)
        
        cd(folders_to_analyze{fold});
        file = pwd;
        %If 0, just do behavior data, no photometry
        start_file = dir(horzcat(file, '/', '*', signal_key_word, '*.csv'));
        if isempty(start_file)
            analyze_photometry_data = 0; 
        else
            analyze_photometry_data = 1; 
        end
        %exponential
        key_words = {bpod_key_word, analog_key_word, signal_key_word, background_key_word};
        AnimalTrackingSaraV7('mov', event_groups, pre_times, post_times, baseline_times, control_fiber, signal_fibers, fiber_names, analyze_photometry_data, key_words)
        %AnimalTrackingSaraV5({16}, {5}, {5}, {5}, 2, 1, fiber_names, analyze_photometry_data, key_words)
    end
    cd(og_file);
    calculate_z_score_w_consistent_baselinesV2();
end