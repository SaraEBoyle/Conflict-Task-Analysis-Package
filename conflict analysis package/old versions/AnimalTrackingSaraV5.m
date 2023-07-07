%% animal trajectory tracking from .avi video
% ---> Inputs: video 
% ---> Outputs: mouse_parameter (column1: Xcenter; column1: Ycenter; column3: Distance; column1: Velocity)
% ---> Sara Boyle sboyle@cshl.edu
%% clear workspace
%% TODO Count up number of seconds spent nose poking
%% TODO- potentially another measure- average point in sound delivery where mouse spends 90% of time on platform
%Idea is that this can be a measure of how well the mouse has learned. If
%less than 18 the mouse has learned
%% TODO- separate platform crosses that are first- ie first platform entry after sound.
%% Find first platform exit after a shock- purple has a long lasting increase in activity to shock. Will she stay on the platform until the signal returns to baseline?
%% First platform exit after sound ends
%% Number of nose pokes
%% Amount of water delivered
%% TODO platform entries not during sound
%% TODO get first water reward
function [] = AnimalTrackingSaraV5(event_groups, pre_times, post_times, baseline_times, control_fiber, signal_fibers, fiber_names, photometry_data, key_words)
%This is a list of all possible events to align to- do not change it
event_labels = {'water', 'light', 'sound', 'laser', 'shock', 'shock received', 'shock avoided', 'sound in shock received', 'sound in shock avoided', 'platform entries', 'platform exits', 'sound platform entries', 'first platform entries', 'first platform exits after shock', 'ITI platform entries', 'first water'};

bpod_key_word = key_words{1};
analog_key_word = key_words{2};
signal_key_word = key_words{3};

record_tracking_video = 0; %if you want to save an example video with the tracking labeled

light_duration = 4;
sound_duration = 20;

%% photometry variables
%This corrects by fitting an exponential
correct_for_bleaching = 0;

%This corrects by using a sliding window
moving_avg = 1;


frame_rate = 30;
fps_photometry = 10;

pre_ITI_length = 150;
seconds_to_extend = 0; 
sound_on = 0;
sound_only = 0;
trigger1 = 1;
plot_single_trials = 0;
plot_heatmap_and_average = 1;
save_graphs = 1;
volume_to_check = 1000;
stim_duration = .2;
sound_offset = 0;
time_after_sound_onset_to_check = 0;
record_sample = 1;
plot_location = 1;
behavior_data = 1;
flip_arena = 1; % just leave this on 1

file = horzcat(pwd, '/');
tem = dir(horzcat(file, '*', signal_key_word, '*.csv')); 
if isempty(tem)
    photometry_data = 0;
end

%% select EthoVision data from excel
% clear nam
if ~exist('nam', 'var') || isempty(nam)
    [file_nm, dir_nm] = uigetfile(fullfile('*.png'));
    if dir_nm==0
        fprintf('no file was selected. Try again.\N');
        
        return;
    end
    nam = [dir_nm, file_nm];  % full name of the data file
    [~, file_nm, file_type] = fileparts(nam);
end

[photometry_start, start_times] = get_start_times(pwd, analog_key_word);

%% Load in data
if photometry_data == 1
    [start_times, photo_data, SessionData, track_data] = load_data(bpod_key_word, signal_key_word, analog_key_word, dir_nm, photometry_data);
end

if size(start_times, 2) > 20
    dif = size(start_times, 2) - 20;
    start_times(1:dif) = [];
    %start_times = start_times - start_times(1);
end

%% Read in video file
if strcmp(file_type, '.mov')
    vidObj = VideoReader([file_nm file_type]); % read data from video
    temp = rgb2gray(readFrame(vidObj));
    %[J,rect2] = imcrop(temp);
    
   
    J = temp;
    rect2 = [0 0 size(temp, 1), size(temp, 2)];

    
    fps = vidObj.FrameRate;
    vidHeight = size(J,1);
    vidWidth = size(J, 2);
    k = 1;
    FrameNumber = vidObj.Duration*vidObj.FrameRate;
    
else
    temp = imread(nam);
    
    J = temp;
    rect2 = [0 0 size(temp, 1), size(temp, 2)];

    fps = 30;
    vidHeight = size(J,1);
    vidWidth = size(J, 2);
    k = 1;
end

file = pwd;
file = horzcat(file, '/');
%video = readframe(vidObj, 1); %reads only the frames specified by index .

close all;
%TODO crop so that you can exclude points not in the box
centroid_file = dir(horzcat(file, '*', 'centroid', '*.csv'));
track_data = dlmread(horzcat(file, centroid_file.name));

if size(track_data, 2) > 4
    centroid_time_scale = track_data(:, 5);
    centroid_time_scale = (centroid_time_scale-photometry_start)/1000;   
else
    centroid_time_scale = 1:size(track_data,1);
    centroid_time_scale = centroid_time_scale./30;
    centroid_time_scale = transpose(centroid_time_scale);
end
new_centroid_time_scale = track_data(:, 5);
new_centroid_time_scale = (new_centroid_time_scale - new_centroid_time_scale(1))/1000;
start_of_bpod_frame = find(centroid_time_scale >= 0);
start_of_bpod_sec = new_centroid_time_scale(start_of_bpod_frame(2));

if isempty(start_of_bpod_frame)
    start_of_bpod_frame = 1;
    disp('issue with tracking');
else
    start_of_bpod_frame = start_of_bpod_frame(1);
end

%start_of_bpod_sec = abs(centroid_time_scale(1));
track_data(1:(start_of_bpod_frame - 1), :) = [];
centroid_time_scale(1:(start_of_bpod_frame - 1)) = [];
%time = new_data(:, 1);
orientation = track_data(:, 3);
length = track_data(:, 4);
track_data(:, 3:end) = [];
track_data = horzcat(track_data, centroid_time_scale);

n = 0.035; % cm/pixel 
mouse_parameter(:,1) = (track_data(:,1) - vidWidth/2)*n;
mouse_parameter(:,2) = (vidHeight/2 - track_data(:,2))*n; 
water_delivered_times = [];
nose_poke_times = [];
light_on_times = [];
light_delivered_ind = [];
laser_on_times = [];
laser_delivered_ind = [];
water_delivered_ind = [];
sound_on_times = [];
sound_delivered_ind = [];
event_end_times = [];
no_event_trial_ind = [];
sound_off_times = [];
shock_times = [];
shock_received_times = [];
shock_avoided_times = [];
shock_avoided_ind = [];
shock_received_ind = [];
sound_in_shock_received = [];
sound_in_shock_avoided = [];
poke_times = [];
first_water_delivered = [];
poke_out_times = [];
poke_time_ind = [];
%Load bpod data
file = dir_nm;
tem = dir(horzcat(file, '*', bpod_key_word, '*.mat')); 
if isempty(tem)
    tem = dir(horzcat(file, '*', bpod_key_word, '*.mat'));
    if isempty(tem)
        tem = dir(horzcat(file, '*behavior*.mat'));
    end
end
if isempty(tem)
    error('set signal_key_word to a phrase that appears in the bpod file')
end

name = tem(1).name;
folder = tem(1).folder;
full_path = horzcat(folder, '/', name);
disp(name);
load(full_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Tracking info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate distance and velocity
for i = 2:size(mouse_parameter,1) % no distance for the first data point
    centerX1 = mouse_parameter(i-1,1);
    centerX2 = mouse_parameter(i,1);
    centerY1 = mouse_parameter(i-1,2);
    centerY2 = mouse_parameter(i,2);
    dis = sqrt((centerX2-centerX1)^2 + (centerY2-centerY1)^2);
    vel = dis/0.2; % cm/s;  sample rate = 5 HZ
    mouse_parameter(i,3) = dis;
    mouse_parameter(i,4) = vel;
end 

%% calculate time spent in each position
new_time_spent = zeros(vidWidth, vidHeight);
sound_time_spent = zeros(vidWidth, vidHeight);
sound_off_time_spent = zeros(vidWidth, vidHeight);
%new_time_spent = zeros(vidWidth, vidHeight);
%sound_time_spent = zeros(vidWidth, vidHeight);
%sound_off_time_spent = zeros(vidWidth, vidHeight);
copy_back = [];
frames_deleted = 0;
first_chunk = 0;
temp;
rect2; % last is first I2, first is last [x, y, width, height] x,y upper left corner of the rectangle
rect2(2); %the first in I2.
rect2(1); %the second in I2.

%new
%d_from_edge_2 = size(temp, 1) - (rect2(2) + rect2(4));
%d_from_edge_1 = size(temp, 2) - (rect2(1) + rect2(3));
d_from_edge_2 = rect2(2);
d_from_edge_1 = rect2(1);
track_data(:, 1) = track_data(:, 1) - d_from_edge_1;
track_data(:, 2) = track_data(:, 2) - d_from_edge_2;

%coords on uncropped image
top_left_corner = [round(rect2(2)), round(rect2(1))];
top_right_corner = [round(rect2(2)), round(rect2(1)) + round(rect2(3))]; %
bottom_left_corner = [round(rect2(2)) + round(rect2(4)), round(rect2(1))];
bottom_right_corner = [round(rect2(2)) + round(rect2(4)), round(rect2(1)) + round(rect2(3))];



%track_data(:, 1) = track_data(:, 1) - top_left_corner(1);
%track_data(:, 2) = track_data(:, 2) - top_left_corner(1);

% pre processing for tracking data

%x_ratio = 1 + ((size(temp, 1) - d_from_edge_2) - size(J, 1))/(size(temp, 1) - d_from_edge_2);
%y_ratio = 1 + ((size(temp, 2) - d_from_edge_1) - size(J, 2))/(size(temp, 2) - d_from_edge_1);
x_ratio = (size(temp, 1) - d_from_edge_2)/size(J, 1);
y_ratio = (size(temp, 2) - d_from_edge_1)/size(J, 2);

for i = 1:size(track_data,1)
    if or(isnan(track_data(i, 1)), isnan(track_data(i, 2)))
        continue
    end
    vert_num = round(track_data(i, 2));
    horz_num = round(track_data(i, 1));
    
    
    
    %I2(bigger = down, smaller = left)

    %I2((abs(top_left_corner(1) - 3)):(top_left_corner(1) + 3), (abs(top_left_corner(2) - 3)):(top_left_corner(2) + 3), 1) = 1;
    %I2((abs(top_right_corner(1) - 3)):(top_right_corner(1) + 3), (abs(top_right_corner(2) - 3)):(top_right_corner(2) + 3), 1) = 1;
    %I2((abs(bottom_right_corner(1) - 3)):(bottom_right_corner(1) + 3), (abs(bottom_right_corner(2) - 3)):(bottom_right_corner(2) + 3), 1) = 1;
    %I2((abs(bottom_left_corner(1) - 3)):(bottom_left_corner(1) + 3), (abs(bottom_left_corner(2) - 3)):(bottom_left_corner(2) + 3), 1) = 1;
    %imshow(I2);
    %This is the spot registered
    %% TODO check ratios
    horz_thresh = round(rect2(3));%* y_ratio);
    vert_thresh = round(rect2(4));% * x_ratio);
    
    if or(horz_num < 0, horz_num > round(rect2(3)))%* y_ratio))
        %if horznum is outside the cropped rectangle
        disp('deleted because y less than 0');
        track_data(i, 2) = NaN;
        track_data(i, 1) = NaN;
        continue
    end
    if or(vert_num < 0, vert_num > round(rect2(4)))% * x_ratio))
        %delete this
        disp('deleted because x less than 0');
        track_data(i, 2) = NaN;
        track_data(i, 1) = NaN;
        continue
    end
    
    if and(vert_num > 3, horz_num > 3)
        %I2 = im2double(J);
        %I2((abs(vert_num - 3)):(vert_num + 3), (abs(horz_num - 3)):(horz_num + 3), 1) = 0;
        %imshow(I2);
    end
    
    %x_num = round(track_data(i, 2)); y_num = round(track_data(i, 1));,
    
    % Find if the mouse is on the platform
    if vert_num < size(J, 1)/2
        %if in top half
        if horz_num > (size(J, 2)/2 + 10)
            %if in right half
            %disp('In platform zone');
            %imshow(I2);
        else
            %disp('Not in platform zone');
        end
    else
        %disp('Not in platform zone');
    end
    %First number is vertical
    %I2(top_left_corner(1) + 20, top_left_corner(2), 1) = 0;
    if vert_num > size(J, 1)/2
        %if in bottom half
        if horz_num  < size(J, 2)/2
            %in left half
            %disp('In reward zone');
            %imshow(I2);
        end
    end
    if or(horz_num + top_left_corner(2) > top_right_corner(2), horz_num < 0)
        %delete this
        disp('weird');
        %track_data(i, 2) = NaN;
        %track_data(i, 1) = NaN;
        continue
    end
    if or(vert_num + top_left_corner(1) > bottom_left_corner(1), vert_num < 0)
        %delete this
        disp('weird');
        %track_data(i, 2) = NaN;
        %track_data(i, 1) = NaN;
        continue
    end
    if or(vert_num >= 4, horz_num >= 4)
        %disp(horzcat('horz_num = ', num2str(horz_num), ', vert_num = ', num2str(vert_num)));
        %I2((abs(vert_num - 3)):(vert_num + 3), (abs(horz_num - 3)):(horz_num + 3), 1) = 0;
        %imshow(I2);
    end
end

%x_ratio = size(temp, 1)/size(J, 1);


%scale x coordinates by top_left_corner(2)
steps = top_left_corner(2)/rect2(3);
%x_cropped_tracking_adjustment = 0:steps:top_left_corner(2);

steps = top_left_corner(1)/rect2(4);
%y_cropped_tracking_adjustment = 0:steps:top_left_corner(1);

%adjusted_x = (1:round(rect2(3)))/top_left_corner(2);
%subtracted = top_left_corner(2) + top_left_corner(1);
%top_right_corner(2) + top_right_corner(2);




%% Preprocess to remove NaN frames
for i = 1:size(track_data,1)
    if and(isnan(track_data(i, 1)), first_chunk == 0)
        %if you're in the first NaN region
        frames_deleted = frames_deleted + 1;
    elseif and(~isnan(track_data(i, 1)), first_chunk == 0)
        %Once you hit a spot where the animal registers for the first time
        first_chunk = 1;
        first_animal_frame = i;
        %disp(horzcat('X point: ', num2str(track_data(i,2)), ' Y point: ', num2str(track_data(i,1)), ' i: ', num2str(i)));
        %new_time_spent(round(track_data(i,1)), round(track_data(i,2))) = new_time_spent(round(track_data(i,1)), round(track_data(i,2))) + 1/fps;
        %CHANGED
        new_time_spent(round(track_data(i,2)), round(track_data(i,1))) = new_time_spent(round(track_data(i,2)), round(track_data(i,1))) + 1/fps;
    else
        if ~isnan(track_data(i, 1))
            if or(round(track_data(i,2)) < 1, round(track_data(i,1)) < 1)
                track_data(i, 1) = NaN;
                track_data(i, 2) = NaN;
                disp('Something is less than 0 here')
            end
        end
        if isnan(track_data(i, 1))
            if i ~= 1
                if isempty(copy_back)
                    track_data(i, 1) = track_data(i - 1, 1);
                    track_data(i, 2) = track_data(i - 1, 2);
                else
                    track_data(i, 1) = 1;
                    track_data(i, 2) = 1;
                    copy_back = horzcat(copy_back, i);
                end
            else
                copy_back = horzcat(copy_back, i);
                track_data(i, 1) = 1;
                track_data(i, 2) = 1;
            end
        elseif ~isempty(copy_back)
            track_data(copy_back, 1) = track_data(i, 1);
            track_data(copy_back, 2) = track_data(i, 2);
            for halp = copy_back
                y_time = round(track_data(halp,2));
                x_time = round(track_data(halp,1));
                new_time_spent(y_time, x_time) = new_time_spent(y_time, x_time) + 1/fps;
                %new_time_spent(x_time, y_time) = new_time_spent(x_time, y_time) + 1/fps;
            end
            copy_back = [];
        end

        if isempty(copy_back)
            %disp(horzcat('X point: ', num2str(track_data(i,2)), ' Y point: ', num2str(track_data(i,1)), ' i: ', num2str(i)));
            %new_time_spent(round(track_data(i,1)), round(track_data(i,2))) = new_time_spent(round(track_data(i,1)), round(track_data(i,2))) + 1/fps;
            %CHANGED
            new_time_spent(round(track_data(i,2)), round(track_data(i,1))) = new_time_spent(round(track_data(i,2)), round(track_data(i,1))) + 1/fps;
            %new_time_spent(round(track_data(i,2)), round(track_data(i,1))) = new_time_spent(round(track_data(i,2)), round(track_data(i,1))) + 1/fps;
        end

        %new_time_spent(round(new_data(i,2)), round(new_data(i,1))) = new_time_spent(round(new_data(i,2)), round(new_data(i,1))) + 1/fps;
    end 
end 


%% Find platform entries
[platform_matrix, sound_data_time_adjusted] = find_sound_platform_position(track_data, J);
platform_entries = [];
platform_exits = [];
prev = platform_matrix(1, 1);
for p = 2:size(platform_matrix, 1)
    if or(p > size(platform_matrix, 1) - 10, p < 11)
        cur = platform_matrix(p, 1);
        prev = cur;
        continue
    end
    cur = platform_matrix(p, 1);
    if cur > prev
        %if mouse crosses to platform
        post = sum(platform_matrix(p:(p+10), 1));
        pre = sum(platform_matrix((p - 11):(p - 1), 1));
        if post >= 8
            %for next 10 values, if most are on platform
            if pre <= 3
                %for previous 10 values, if most are not on platform
                platform_entries = horzcat(platform_entries, platform_matrix(p, 2));
                prev = cur;
            end
        end
    elseif cur < prev
        post = sum(platform_matrix(p:(p+10), 1));
        pre = sum(platform_matrix((p - 11):(p - 1), 1));
        %if mouse crosses from platform
        if post <= 3
            %if most are not on platform
            if pre >= 8
                %for previous 10 values, if most are on platform
                platform_exits = horzcat(platform_exits, platform_matrix(p, 2));
                prev = cur;
            end
        end
    end
    prev = cur;
    
end






%Find time of events

for q = 1:min(SessionData.nTrials, size(start_times, 2))
    trial_start = start_times(q);
    poke_port = 'Port2In';
    reward_state = 'Reward';
    %Record nose poke times
    if isfield(SessionData.RawEvents.Trial{1, q}.Events, poke_port)
        if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.Events.', poke_port, '(1)')))
            poke_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.Events.', poke_port));
            if ~isnan(poke_time)        
                nose_poke_times = horzcat(nose_poke_times, poke_time + trial_start);
            end
        end
    end
    
    
    
    %Record water delivery times
    if isfield(SessionData.RawEvents.Trial{1, q}.States, reward_state)
        if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', reward_state, '(1)')))
            reward_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', reward_state));
            reward_time = reward_time(:,1);
            if ~isnan(reward_time)        
                water_delivered_times = vertcat(water_delivered_times, reward_time + trial_start);
                water_delivered_ind = horzcat(water_delivered_ind, q); 
                first_water_delivered = vertcat(first_water_delivered, reward_time(1) + trial_start);
            end
        end
    end

    %Record light delivery times
    if isfield(SessionData.RawEvents.Trial{1, q}.States, 'LightCueOn')
        if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', 'LightCueOn', '(1)')))
            light_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', 'LightCueOn'));
            light_on_time = light_time(:,1);
            light_off_time = light_time(:,2);

            if ~isnan(light_time)        
                light_on_times = vertcat(light_on_times, light_on_time + trial_start);
                light_delivered_ind = horzcat(light_delivered_ind, q);               
            end
        end
    end
    
    %Record poke in
    if isfield(SessionData.RawEvents.Trial{1, q}.Events, 'Port2In')
        if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.Events.', 'Port2In', '(1)')))
            poke_in_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.Events.', 'Port2In'));
            poke_out_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.Events.', 'Port2Out'));
            if ~isnan(poke_in_time)  
                
                if size(poke_in_time, 2) > size(poke_out_time, 2)
                    if poke_in_time(1) > poke_out_time(1)
                        poke_in_time(1) = [];
                    else
                        poke_in_time(end) = [];
                    end
                elseif size(poke_in_time, 2) < size(poke_out_time, 2)
                    if poke_in_time(1) > poke_out_time(1)
                        poke_out_time(1) = [];
                    elseif poke_in_time(end) > poke_out_time(end)
                        poke_out_time(end) = [];
                    else
                        poke_out_time(end) = [];
                    end
                end
                
                pokes =  vertcat(poke_in_time, poke_out_time);
                
                poke_times = horzcat(poke_times, pokes + trial_start);
                %poke_times = vertcat(poke_out_times, poke_out_time + trial_start);              
                poke_time_ind = horzcat(poke_time_ind, q);               
            end
            
            

            %if ~isnan(poke_out_time)        
            %    poke_out_times = vertcat(poke_out_times, poke_out_time + trial_start);              
            %end
        end
    end
    
    %Record laser delivery times
    if isfield(SessionData.RawEvents.Trial{1, q}.States, 'Laser')
        if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', 'Laser', '(1)')))
            laser_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', 'Laser'));
            laser_on_time = laser_time(:,1);

            if ~isnan(laser_time)        
                laser_on_times = vertcat(laser_on_times, laser_on_time + trial_start);
                laser_delivered_ind = horzcat(laser_delivered_ind, q);               
            end
        end
    end

    %Record sound delivered times
    %Find sound event-
    if isfield(SessionData.RawEvents.Trial{1, q}.States, 'SoundOn')
        sound_state = 'SoundOn';
    elseif isfield(SessionData.RawEvents.Trial{1, q}.States, 'SoundCueOn2')
        sound_state = 'SoundCueOn2';
    elseif isfield(SessionData.RawEvents.Trial{1, q}.States, 'SoundCueOn')
        sound_state = 'SoundCueOn';
    end
    
    if isfield(SessionData.RawEvents.Trial{1, q}.States, sound_state)
        if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, ', num2str(q), '}.States.', sound_state, '(1)')))
            stimulus_time = eval(horzcat('SessionData.RawEvents.Trial{1, ', num2str(q), '}.States.', sound_state));
            sound_on_time = stimulus_time(1);
            sound_off_time = stimulus_time(2);
            sound_on_times = horzcat(sound_on_times, sound_on_time + trial_start);
            sound_off_times = horzcat(sound_off_times, sound_off_time + trial_start);
            sound_delivered_ind = horzcat(sound_delivered_ind, q);
            
            %Sort out shock
            shock_times = horzcat(shock_times, sound_on_time + 18 + trial_start);
            shock_on = sound_on_time + 18 + trial_start;
            first_shock_pos = find(track_data(:,3) >= shock_on);
            first_shock_pos = track_data(first_shock_pos(1), :);
            vert_num = round(first_shock_pos(1, 2));
            horz_num = round(first_shock_pos(1, 1));
            
            if ~flip_arena
                % Find if the mouse is on the platform (top right)
                if vert_num < size(J, 1)/2
                    %if in top half
                    if horz_num > (size(J, 2)/2 + 10)
                        %if in right half
                        %disp('In platform zone');
                        %imshow(I2);
                        disp('Shock avoided');
                        shock_avoided_times = horzcat(shock_avoided_times, sound_on_time + 18 + trial_start);
                        sound_in_shock_avoided = horzcat(sound_in_shock_avoided, sound_on_time + trial_start);
                        shock_avoided_ind = horzcat(shock_avoided_ind, q);
                    else
                        disp('Shock received');
                        %disp('Not in platform zone');
                        shock_received_times = horzcat(shock_received_times, sound_on_time + 18 + trial_start);
                        sound_in_shock_received = horzcat(sound_in_shock_received, sound_on_time + trial_start);
                        shock_received_ind = horzcat(shock_received_ind, q);
                    end
                else
                    %disp('Not in platform zone');
                    disp('Shock received');
                    shock_received_times = horzcat(shock_received_times, sound_on_time + 18 + trial_start);
                    sound_in_shock_received = horzcat(sound_in_shock_received, sound_on_time + trial_start);
                    shock_received_ind = horzcat(shock_received_ind, q);
                end
            else
                %if platform is bottom left (flip arena)
                %I2(top_left_corner(1) + 20, top_left_corner(2), 1) = 0;
                if vert_num > size(J, 1)/2
                    %if in bottom half
                    if horz_num  < size(J, 2)/2
                        %in left half
                        disp('Shock avoided');
                        %imshow(I2);
                        shock_avoided_times = horzcat(shock_avoided_times, sound_on_time + 18 + trial_start);
                        sound_in_shock_avoided = horzcat(sound_in_shock_avoided, sound_on_time + trial_start);
                        shock_avoided_ind = horzcat(shock_avoided_ind, q);
                    else
                        disp('Shock received');
                        shock_received_times = horzcat(shock_received_times, sound_on_time + 18 + trial_start);
                        sound_in_shock_received = horzcat(sound_in_shock_received, sound_on_time + trial_start);
                        shock_received_ind = horzcat(shock_received_ind, q);
                    end
                else
                    disp('Shock received');
                    shock_received_times = horzcat(shock_received_times, sound_on_time + 18 + trial_start);
                    sound_in_shock_received = horzcat(sound_in_shock_received, sound_on_time + trial_start);
                    shock_received_ind = horzcat(shock_received_ind, q);
                end
            end
            
        else
            sound_on_times = horzcat(sound_on_times, 0); %Record trials with no sound as the start of the trial
        end
    end
    
end

%% Find first platform entries during sound deliveries
platform_entries_during_sound = [];
first_platform_entries_during_sound = [];
platform_entries_during_ITI = [];
first = 0;
for q = 1:size(sound_on_times, 2)
    for m = platform_entries
        if m < sound_on_times(q)
            %if entry is before sound onset skip
            continue
        elseif m <= sound_on_times(q) + 18
            %If entry is after sound onset and before shock onset
            platform_entries_during_sound = horzcat(platform_entries_during_sound, m);
            if first == 0
                first_platform_entries_during_sound = horzcat(first_platform_entries_during_sound, m);
                first = m;
            end
        elseif q + 1 > size(sound_on_times, 2)
            break
        elseif and(m < sound_on_times(q + 1), m > sound_on_times(q) + 20)
            %if entry is before next sound and after this sound
            platform_entries_during_ITI = horzcat(platform_entries_during_ITI, m);
        
        elseif m >= sound_on_times(q + 1)
            %if entry is after start of next sound move on to next entry
            break
        end
    end
    first = 0;
end
%% %% Find first platform exit after shock deliveries
platform_exits_after_shock = [];
first_platform_exits_after_shock = [];
first = 0;
for q = shock_received_times
    for m = platform_exits
        if m < q
            %if entry is before sound onset skip
            continue
        elseif m <= q + 90
            %If entry is after sound onset and before shock onset
            platform_exits_after_shock = horzcat(platform_exits_after_shock, m);
            if first == 0
                first_platform_exits_after_shock = horzcat(first_platform_exits_after_shock, m);
                first = m;
            end
        elseif m > q + 90
            %if entry is after end of sound check next sound time
            break
        end
    end
    first = 0;
end

%photometry data isn't pulled in until here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

event_times = {};
event_times{1} = water_delivered_times;
event_times{2} = light_on_times;
event_times{3} = transpose(sound_on_times);
event_times{4} = laser_on_times;
event_times{5} = transpose(shock_times);
event_times{6} = transpose(shock_received_times);
event_times{7} = transpose(shock_avoided_times);
event_times{8} = transpose(sound_in_shock_received);
event_times{9} = transpose(sound_in_shock_avoided);
event_times{10} = transpose(platform_entries);
event_times{11} = transpose(platform_exits);
event_times{12} = transpose(platform_entries_during_sound);
event_times{13} = transpose(first_platform_entries_during_sound);
event_times{14} = transpose(first_platform_exits_after_shock);
event_times{15} = transpose(platform_entries_during_ITI);
event_times{16} = first_water_delivered;

save(horzcat(file, 'event_times.mat'), 'water_delivered_times', 'light_on_times', ...
    'sound_on_times', 'laser_on_times', 'shock_times', 'shock_received_times', ...
    'shock_avoided_times', 'sound_in_shock_received', 'sound_in_shock_avoided', ...
    'platform_entries', 'platform_exits', 'platform_entries_during_sound', ...
    'first_platform_entries_during_sound', 'first_platform_exits_after_shock', ...
    'platform_entries_during_ITI', 'first_water_delivered');

raw_signal_dFFs = {};
raw_control_dFFs = {};
event_group_ind = 1;

for ev_ind = 1:size(event_groups, 2)
    
    events_to_analyze = event_groups{ev_ind};
    pre_time = pre_times{ev_ind};
    post_time = post_times{ev_ind};
    baseline_time = baseline_times{ev_ind};
    for fiber = signal_fibers
        seconds_per_trial = pre_time + post_time;
        fiber_name = fiber_names{fiber};
        if photometry_data == 1
            %%So it's the event time minus the first time value in data to get the
            %%timestamp in the video

            for signal = 1:2
                [~, photo_data, ~, ~] = load_data(bpod_key_word, signal_key_word, analog_key_word, dir_nm, photometry_data);
                if signal == 1
                    fiber_name = horzcat(fiber_names{fiber}, ' signal');
                    %to get signal data
                    %Trim your data if needed
                    [photo_data, auto_flo] = trim_data(photo_data, fiber, 1);

                    %De-interleave and subtract autofluorescence
                    [gcamp, iso] = de_interleave(1, 0, file, fiber_name, photo_data, fiber, auto_flo, 1, signal);

                    %correct bleaching
                    [gcamp, iso] = correct_bleaching(correct_for_bleaching, gcamp, iso, moving_avg, frame_rate, 1, 0);

                elseif signal == 2
                    %to get control channel
                    fiber_name = horzcat(fiber_names{fiber}, ' control');
                    if fiber == control_fiber
                        %if it's the same channel as signal
                        gcamp = iso; %just swap the channels
                    else
                        %otherwise you have to re analyze
                        fiber = control_fiber;
                        %Trim your data if needed

                        %signal_file = dir(horzcat(dir_nm, '*', signal_key_word, '*.csv'));
                        %photo_data = dlmread(horzcat(dir_nm, signal_file.name));
                        [IAmLazy, photo_data, alsotrash, trash] = load_data(bpod_key_word, signal_key_word, analog_key_word, dir_nm, photometry_data);

                        [photo_data, auto_flo] = trim_data(photo_data, fiber, 1);

                        %De-interleave and subtract autofluorescence
                        [gcamp, iso] = de_interleave(1, 0, file, fiber_name, photo_data, fiber, auto_flo, 1, signal);

                        %correct bleaching
                        [gcamp, iso] = correct_bleaching(correct_for_bleaching, gcamp, iso, moving_avg, frame_rate, 1, 0);

                    end
                    
                end
                if and(event_group_ind == 1, signal == 1)
                    raw_signal_dFFs{fiber} = gcamp;
                elseif and(event_group_ind == 1, signal == 2)
                    raw_control_dFFs{fiber} = gcamp;
                end

                for event_type_ind = events_to_analyze
                    if isempty(event_times{event_type_ind})
                        continue
                    else
                        alignment = event_labels{event_type_ind};     
                    end
                    if ~exist(horzcat(file, fiber_name), 'dir')
                        mkdir(horzcat(file, fiber_name));
                        mkdir(horzcat(file, fiber_name, '/', alignment));
                    elseif ~exist(horzcat(file, fiber_name, '/', alignment))
                        mkdir(horzcat(file, fiber_name, '/', alignment));
                    end
                    filesub = horzcat(file, fiber_name, '/', alignment, '/');
                    cd(filesub);


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
                    event_frame = {};          %Record frame air puff happens                              
                    event_delivered = [];

                    %Just for no US trials
                    no_event_trials = {};
                    no_event_i = 1;
                    no_event_pre_ITIs = {};
                    no_event_post_ITIs = {};
                    sound_delivered_times_temp = event_times{3};
                    if strcmp(alignment, 'sound')
                        times_delivered = event_times{3};
                        event_times_temp = event_times{3};
                        event_trial_ind = sound_delivered_ind;
                    elseif strcmp(alignment, 'water')
                        times_delivered = event_times{1};
                        event_times_temp = event_times{1};
                        event_trial_ind = water_delivered_ind;
                    elseif strcmp(alignment, 'laser')
                        times_delivered = event_times{4};
                        event_times_temp = event_times{4};
                        event_trial_ind = laser_delivered_ind;
                    elseif strcmp(alignment, 'light')
                        times_delivered = event_times{2};
                        event_times_temp = event_times{2};
                        event_trial_ind = light_delivered_ind;
                    elseif strcmp(alignment, 'shock')
                        times_delivered = event_times{5};
                        event_times_temp = event_times{5};
                        event_trial_ind = sound_delivered_ind;
                    elseif strcmp(alignment, 'shock received')
                        times_delivered = event_times{6};
                        event_times_temp = event_times{6};
                        event_trial_ind = shock_received_ind;
                    elseif strcmp(alignment, 'shock avoided')
                        times_delivered = event_times{7};
                        event_times_temp = event_times{7};
                        event_trial_ind = shock_avoided_ind;
                    elseif strcmp(alignment, 'sound in shock received')
                        times_delivered = event_times{8};
                        event_times_temp = event_times{8};
                        event_trial_ind = shock_received_ind;
                    elseif strcmp(alignment, 'sound in shock avoided')
                        times_delivered = event_times{9};
                        event_times_temp = event_times{9};
                        event_trial_ind = shock_avoided_ind;
                    elseif strcmp(alignment, 'platform entries')
                        times_delivered = event_times{10};
                        event_times_temp = event_times{10};
                        event_trial_ind = 1:20;
                    elseif strcmp(alignment, 'platform exits')
                        times_delivered = event_times{11};
                        event_times_temp = event_times{11};
                        event_trial_ind = 1:20;
                    elseif strcmp(alignment, 'sound platform entries')
                        times_delivered = event_times{12};
                        event_times_temp = event_times{12};
                        event_trial_ind = 1:20; 
                    elseif strcmp(alignment, 'first platform entries')
                        times_delivered = event_times{13};
                        event_times_temp = event_times{13};
                        event_trial_ind = 1:20; 
                    elseif strcmp(alignment, 'first platform exits after shock')
                        times_delivered = event_times{14};
                        event_times_temp = event_times{14};
                        event_trial_ind = 1:20; 
                    elseif strcmp(alignment, 'ITI platform entries')
                        times_delivered = event_times{15};
                        event_times_temp = event_times{15};
                        event_trial_ind = 1:20;   
                    elseif strcmp(alignment, 'first water')
                        times_delivered = event_times{16};
                        event_times_temp = event_times{16};
                        event_trial_ind = water_delivered_ind;
                    end

                    %start_times = alt_start_times; 
                    start_times = times_delivered - pre_time;
                    end_times = times_delivered + post_time;
                    evs = 1;
                    %TODO INVESTIGATE
                    gcamp_clone = gcamp;
                    for y = 1:size(start_times, 1)
                        %Segment the trials based on start and end times
                        %pre_ITI_length = 5;
                        post_ITI_length = 150;   
                        %if evs > length(event_times)

                        starty = start_times(y);
                        endy = end_times(y);
                        indices1 = find(gcamp(:,1) >= starty);
                        if 1
                            indices2 = find(iso(:,1) >= starty);
                            trial2= iso(indices2,:);
                        end
                        indices3 = find(gcamp_clone(:,1) >= starty);
                        trial1=gcamp(indices1,:);
                        trial3 = gcamp_clone(indices3, :);

                        indices1 = find(trial1(:,1) <= endy);
                        if 1
                            indices2 = find(trial2(:,1) <= endy);
                            trial2=trial2(indices2,:);
                        end
                        indices3 = find(trial3(:,1) <= endy);
                        trial1=trial1(indices1,:);
                        if size(trial1, 1) < fps_photometry*seconds_per_trial - 2
                            event_times_temp(1) = [];
                            disp('Event too soon. Had to skip first event')
                            continue
                        end
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
                            else
                                if size(trial1, 1) < ((pre_time + post_time) * 9.5)
                                    event_times_temp(1) = [];
                                    disp('recording is off. Had to skip first event');
                                    continue
                                end
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
                            US_delivered = horzcat(US_delivered, y);         %record trial # with US
                            trials{ind} = trial1(:, [1 2]);                  %Record trial data in total list
                            if trigger1
                                iso_trials{ind} = trial2(:, [1 2]);
                            end
                            OG_trials{ind} = trial3(:, [1 2]);
                            pre_ITIs{ind} = pre_ITI;                         %Record pre ITI data
                            post_ITIs{ind} = post_ITI;                       %Record post ITI data

                            if ismember(y, US_delivered)
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
                            if ind == 0
                                ind = 1;
                            end
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
                    %        checked_trial_types = horzcat(checked_trial_types, trial_types(event_trial_ind(y)));
                            pre_ITIs{ind} = pre_ITI;
                            post_ITIs{ind} = post_ITI;
                            ind = ind + 1;
                        end
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
                        if ~exist('Alignment Data/', 'dir')
                            mkdir(horzcat(file, 'Alignment Data/'));
                        end
                        align_to_this = alignment;
                        if ~isempty(event_trial_ind)
                            good_pre_ITIs = pre_ITIs;
                            good_post_ITIs = post_ITIs;

                            [US_trials, US_target_frame, true_seconds_per_trial, sec_offsets, true_sec_offsets] = sort_and_align_trials_V7(event_trials, good_pre_ITIs, good_post_ITIs, ...
                                seconds_per_trial, US_frame, 10, seconds_to_extend);

                            frame_offsets = sec_offsets * 10;
                            trial_start_times = start_times;
                            trial_order = 1:size(start_times, 1);
                            save(horzcat(file, 'Alignment Data/', alignment, ' alignment data.mat'), 'sec_offsets', 'frame_offsets', 'event_trial_ind', 'no_event_trial_ind', 'true_sec_offsets', 'event_times', 'trial_start_times', 'trial_order');
                        end
                        no_US_trials = [];
                        if ~isempty(event_trial_ind)
                            US_trials = US_trials * 100;
                            baselines = [];
                            for r = 1:size(US_trials, 1)
                                cur = US_trials(r, :);
                                baseline = cur(1:((baseline_time - 1)*fps_photometry - 5));
                                baselines(r, :) = baseline;
                            end
                            plot_average('010', fiber_name, US_trials, true_seconds_per_trial, US_target_frame, 0, 0, alignment, .2);
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
                                    savefig(horzcat(alignment, '_average.fig'));
                                end
                            end
                        else
                            disp(horzcat('Found no ', alignment, ' trials. Exited.'))
                            cd(file);

                        end

                        %% Plot trial by trial responses


                        baseline_no_US_trials = [];
                        las_baseline_no_US_trials = [];
                        no_las_baseline_no_US_trials = [];
                        z_no_US_trials = [];
                        %end
                        if ~isempty(event_trial_ind)            
                            plot_heatmap(US_trials, fiber_name, true_seconds_per_trial, US_target_frame, 0, alignment, pre_time);
                            title(horzcat(fiber_name, ': ', alignment, ' trials'));

                            if strcmp(alignment, 'left')
                                savefig('left_delivery_heatmap.fig');
                            elseif strcmp(alignment, 'right')
                                savefig('right_delivery_heatmap.fig');
                            else
                                savefig(horzcat(alignment, '_heatmap.fig'));
                            end              
                        end

                        increments = size(US_trials, 2);
                        if increments == 0
                            increments = size(no_US_trials, 2);
                        end
                        x_axis = ((1:increments)/increments)*true_seconds_per_trial;
                        if volume_to_check == 1000
                            laser_time = x_axis(US_target_frame);
                            x_axis = x_axis - laser_time;
                            if ~isempty(event_trial_ind)
                                if strcmp(alignment, 'left')
                                    left_trials = US_trials;
                                    save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                elseif strcmp(alignment, 'right')
                                    right_trials = US_trials;
                                    save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                else
                                    save('dFF_data.mat', 'US_trials', 'US_target_frame', 'x_axis', 'baselines', 'auto_flo', 'pre_time', 'post_time');
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
                                plot_average(010, fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
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
                                    plot_average(010, fiber_name, z_no_US_trials, true_seconds_per_trial, CS_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
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

                                plot_average('010', fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, 0, 0, alignment, stim_duration);
                                title(horzcat(fiber_name, ': ', alignment, ' trials'));
                                ylabel('z-score(dF/F)', 'FontSize', 20);
                                if strcmp(alignment, 'left')
                                    savefig('left_z_score_average.fig');
                                    save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                    z_left_trials = z_US_trials;
                                    save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                                elseif strcmp(alignment, 'right')
                                    savefig('right_z_score_average.fig');
                                    save('dFF_data.mat', 'right_trials', 'US_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                    z_right_trials = z_US_trials;
                                    save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                                else
                                    %save('dFF_data.mat', 'right_trials', 'US_target_frame', 'x_axis', 'baselines');
                                    save(horzcat(filesub, 'z_scores.mat'), 'z_US_trials', 'x_axis', 'US_target_frame');
                                    z_right_trials = z_US_trials;
                                    savefig(horzcat(filesub, 'US_z_score_average.fig'));
                                end
                                if and(~isempty(no_event_trial_ind), ~isempty(no_US_trials))
                                    plot_average(010, fiber_name, z_no_US_trials, true_seconds_per_trial, CS_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                                    title(horzcat(fiber_name, ': ', alignment, ' trials'));
                                    ylabel('z-score(dF/F)', 'FontSize', 20);
                                    if strcmp(alignment, 'left') 
                                        save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                        z_right_trials = z_no_US_trials;
                                        z_left_trials = z_US_trials;
                                        save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                                    elseif strcmp(alignment, 'right')
                                        save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
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
                                save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                            elseif strcmp(alignment, 'left')
                                save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                            else
                                save('dFF_data.mat', 'US_trials', 'no_US_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'auto_flo', 'pre_time', 'post_time');
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
                                plot_average(010, fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
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
            end
        %random_message(2);
        end
    end
    event_group_ind = event_group_ind + 1;
end

save(horzcat(file, 'raw_data.mat'), 'raw_signal_dFFs', 'raw_control_dFFs');
sound_on_nose_pokes = 0;
sound_off_nose_pokes = 0;
sound_on_times_copy = sound_on_times;
sound_on_poke_time = 0;
sound_off_poke_time = 0;

for s = 1:size(poke_times, 2)
    if and(nose_poke_times(1, s) > sound_on_times_copy(1), nose_poke_times(1, s) < sound_on_times_copy(1) + sound_duration - sound_offset)
        %If after sound but before end of sound
        sound_on_nose_pokes = sound_on_nose_pokes + 1;
        sound_on_poke_time = sound_on_poke_time + abs(poke_times(2, s) - poke_times(1, s));
    else
        %If not in sound window
        if (sound_on_times_copy(1) + sound_duration - sound_offset) < nose_poke_times(1, s)
            %If after end of sound
            sound_off_nose_pokes = sound_off_nose_pokes + 1;
            sound_off_poke_time = sound_off_poke_time + abs(poke_times(2, s) - poke_times(1, s));
            sound_on_times_copy(1) = [];
            if isempty(sound_on_times_copy)
                break
            end
        elseif (sound_on_times_copy(1) - sound_offset) > nose_poke_times(1, s)
            %if before start of sound
            sound_off_nose_pokes = sound_off_nose_pokes + 1;
            sound_off_poke_time = sound_off_poke_time + abs(poke_times(2, s) - poke_times(1, s));
        end
    
    end
    if (sound_on_times_copy(1) + sound_duration - sound_offset) < nose_poke_times(1, s)
        %If after end of sound
        sound_on_times_copy(1) = [];
        if isempty(sound_on_times_copy)
            break
        end
    end
end

session_start = SessionData.TrialStartTimestamp(1);
session_end = SessionData.TrialEndTimestamp(end);
session_length = session_end - session_start;
total_sound_time = SessionData.nTrials * sound_duration;
total_sound_off_time = session_length - total_sound_time;

if ~exist('Behavior')
    mkdir('Behavior');
end

%plot number of pokes
X = categorical({'Total', 'Sound Off', 'Sound On'});
Y = [(sound_off_nose_pokes + sound_on_nose_pokes)/session_length, sound_off_nose_pokes/total_sound_off_time, sound_on_nose_pokes/total_sound_time];
bar(X, Y);
title('Number of nose pokes per second');
savefig('Behavior/number_of_nose_pokes.fig')

%Plot time spent poking
X = categorical({'Total', 'Sound Off', 'Sound On'});
X = reordercats(X, {'Total', 'Sound Off', 'Sound On'});
all_poke_time = sound_on_poke_time + sound_off_poke_time;
c_all_poke_time = all_poke_time/session_length;
c_sound_off_poke_time = sound_off_poke_time/total_sound_off_time;
c_sound_on_poke_time = sound_on_poke_time/total_sound_time;
Y = [all_poke_time, sound_off_poke_time, sound_on_poke_time];
bar(X, Y);
title('Nose poke duration');
savefig('Behavior/nose_poke_duration.fig')
save('Behavior/nose_poke_duration.mat', 'all_poke_time', 'sound_on_poke_time', 'sound_off_poke_time');

%Plot %time spent poking
Y = [c_all_poke_time, c_sound_off_poke_time, c_sound_on_poke_time];
bar(X, Y);
title('Nose poke duration');
ylabel('% time spent')
savefig('Behavior/nose_poke_duration_percent.fig')
save('Behavior/nose_poke_duration_percent.mat', 'c_all_poke_time', 'c_sound_on_poke_time', 'c_sound_off_poke_time');


total_shocks_received = size(shock_received_times, 2);
total_shocks_avoided = size(shock_avoided_times, 2);
total_water_delivered = size(water_delivered_times, 1);
save('Behavior/session_performance', 'total_shocks_received', 'total_shocks_avoided', 'total_water_delivered', 'sound_on_poke_time', 'sound_off_poke_time')



if plot_location

    %% get location during the sound delivery
    sound_on_times_copy = sound_on_times;

    %Delete everything except the sound on period
    sound_on_data = [];
    size_sound = 1;
    for k = 1:size(track_data,1)
        if and(track_data(k, 3) > sound_on_times_copy(1) + time_after_sound_onset_to_check, track_data(k, 3) < sound_on_times_copy(1) + sound_duration - sound_offset)
            %If after sound but before end of sound
            sound_on_data(size_sound, :) = track_data(k, :);
            size_sound = size_sound + 1;
        end
        if (sound_on_times_copy(1) + sound_duration - sound_offset) < track_data(k, 3)
            %If after end of sound
            sound_on_times_copy(1) = [];
            if isempty(sound_on_times_copy)
                break
            end
        end
    end

    % Collect the locations during sound on
    copy_back = [];
    for i = 1:size(sound_on_data,1)
        %TODO segment by the events. Light on, sound on, total
        if isnan(sound_on_data(i, 1))
            if i ~= 1
                if isempty(copy_back)
                    sound_on_data(i, 1) = sound_on_data(i - 1, 1);
                    sound_on_data(i, 2) = sound_on_data(i - 1, 2);
                else
                    sound_on_data(i, 1) = 1;
                    sound_on_data(i, 2) = 1;
                    copy_back = horzcat(copy_back, i);
                end
            else
                copy_back = horzcat(copy_back, i);
                sound_on_data(i, 1) = 1;
                sound_on_data(i, 2) = 1;
            end
        elseif ~isempty(copy_back)
            sound_on_data(copy_back, 1) = sound_on_data(i, 1);
            sound_on_data(copy_back, 2) = sound_on_data(i, 2);
            for halp = copy_back
                y_time = round(sound_on_data(halp,2));
                x_time = round(sound_on_data(halp,1));
                %sound_time_spent(x_time, y_time) = sound_time_spent(x_time, y_time) + 1/fps;
                sound_time_spent(y_time, x_time) = sound_time_spent(y_time, x_time) + 1/fps;
            end
            copy_back = [];
        end

        if isempty(copy_back)
            %sound_time_spent(round(sound_on_data(i,1)), round(sound_on_data(i,2))) = sound_time_spent(round(sound_on_data(i,1)), round(sound_on_data(i,2))) + 1/fps;
            sound_time_spent(round(sound_on_data(i,2)), round(sound_on_data(i,1))) = sound_time_spent(round(sound_on_data(i,2)), round(sound_on_data(i,1))) + 1/fps;
        end
    end 

    %% get position while sound is off
    sound_on_times_copy = sound_on_times;

    %Delete everything except the sound on period
    sound_off_data = [];
    size_sound = 1;
    trial_num = 1;
    for k = 1:size(track_data,1)
        if trial_num > size(start_times, 1)
            break
        end
        if isfield(SessionData.RawEvents.Trial{1, trial_num}.States, 'ITI')
            ITI_start = SessionData.RawEvents.Trial{1, trial_num}.States.ITI(1) + start_times(trial_num);
            ITI_end = SessionData.RawEvents.Trial{1, trial_num}.States.ITI(2) + start_times(trial_num);
        else
            ITI_start = SessionData.RawEvents.Trial{1, trial_num}.Events.GlobalTimer3_End + start_times(trial_num);
            ITI_end = SessionData.RawEvents.Trial{1, trial_num}.Events.GlobalTimer3_End + start_times(trial_num);
        end
        
        if ~and(track_data(k, 3) > sound_on_times_copy(1), track_data(k, 3) < (sound_on_times_copy(1) + sound_duration - sound_offset))
            %If NOT after sound start and before sound end
            if ~and(track_data(k, 3) > ITI_start, track_data(k, 3)<  ITI_end)
                %If NOT within ITI
                sound_off_data(size_sound, :) = track_data(k, :);
                size_sound = size_sound + 1;
            end
        end
        %if (sound_on_times_copy(1) + sound_duration - sound_offset) < track_data(k, 3)
        if ITI_end < track_data(k, 3)
            %if after the end of the trial
            sound_on_times_copy(1) = [];
            trial_num = trial_num + 1;
            if isempty(sound_on_times_copy)
                break
            end
        end
    end

    % Collect the locations during sound on
    copy_back = [];
    for i = 1:size(sound_off_data,1)
        %TODO segment by the events. Light on, sound on, total
        if isnan(sound_off_data(i, 1))
            if i ~= 1
                if isempty(copy_back)
                    sound_off_data(i, 1) = sound_off_data(i - 1, 1);
                    sound_off_data(i, 2) = sound_off_data(i - 1, 2);
                else
                    sound_off_data(i, 1) = 1;
                    sound_off_data(i, 2) = 1;
                    copy_back = horzcat(copy_back, i);
                end
            else
                copy_back = horzcat(copy_back, i);
                sound_off_data(i, 1) = 1;
                sound_off_data(i, 2) = 1;
            end
        elseif ~isempty(copy_back)
            sound_off_data(copy_back, 1) = sound_off_data(i, 1);
            sound_off_data(copy_back, 2) = sound_off_data(i, 2);
            for halp = copy_back
                y_time = round(sound_off_data(halp,2));
                x_time = round(sound_off_data(halp,1));
                %sound_off_time_spent(x_time, y_time) = sound_off_time_spent(x_time, y_time) + 1/fps;
                sound_off_time_spent(y_time, x_time) = sound_off_time_spent(y_time, x_time) + 1/fps;
            end
            copy_back = [];
        end

        if isempty(copy_back)
            %sound_off_time_spent(round(sound_off_data(i,1)), round(sound_off_data(i,2))) = sound_off_time_spent(round(sound_off_data(i,1)), round(sound_off_data(i,2))) + 1/fps;
            sound_off_time_spent(round(sound_off_data(i,2)), round(sound_off_data(i,1))) = sound_off_time_spent(round(sound_off_data(i,2)), round(sound_off_data(i,1))) + 1/fps;
        end
    end 


    %% Smooth the tracking data
    smoothed_data_x = track_data(:,1);
    smoothed_data_y = track_data(:,2);
    %smoothed_data_x(1:(first_animal_frame - 1), :) = [];
    %smoothed_data_y(1:(first_animal_frame - 1), :) = [];
    %smoothed_track_data = track_data;
    %%{
    smoothed_data_x = smooth(smoothed_data_x);
    smoothed_data_y = smooth(smoothed_data_y);
    smoothed_track_data = horzcat(smoothed_data_x, smoothed_data_y);
    %%}
    %smoothed_new_data(1:(first_animal_frame - 1), :) = [];
    nose_points = track_data;
    [x,y,z] = sph2cart(orientation,zeros(1, size(orientation, 2)),1);
    %nose_points(:, 1) = smooth(30.*x);
    %nose_points(:, 2) = smooth(30.*y);
    %% TODO divide into quadrants and quantify time on platform vs reward area
    %I switched the order so now i need to change rows and columns
    %all time
    rows = round(size(new_time_spent, 2)/2);
    columns = round(size(new_time_spent, 1)/2);
    top_left = new_time_spent(1:rows, 1:columns);
    bottom_left = new_time_spent((rows + 1):end, 1:columns);
    top_right = new_time_spent(1:rows, (columns + 1):end);
    bottom_right = new_time_spent((rows + 1):end, (columns + 1):end);
    
    I2 = im2double(J);
    I2((abs(rows - 3)):(rows + 3), (abs(columns - 3)):(columns + 3), 1) = 0;
    I2((abs(4 - 3)):(4 + 3), (abs(4 - 3)):(4 + 3), 1) = 0;
    I2((abs(rows - 3)):(rows + 3), (abs(columns*2 - 3)):(columns*2), 1) = 0;
    imshow(I2);
    
    if ~flip_arena
        top_left_time = sum(sum(top_left));
        platform_zone_time = sum(sum(top_right));
        reward_zone_time = sum(sum(bottom_left));
        bottom_right_time = sum(sum(bottom_right));
    else
        top_left_time = sum(sum(top_left));
        reward_zone_time = sum(sum(top_right));
        platform_zone_time = sum(sum(bottom_left));
        bottom_right_time = sum(sum(bottom_right));
        %changed
        %top_left_time = sum(sum(bottom_right));
        %platform_zone_time = sum(sum(bottom_left));
        %reward_zone_time = sum(sum(top_right));
        %bottom_right_time = sum(sum(top_left));
    end
    %all_percent_reward_zone = reward_zone_time/(top_left_time + platform_zone_time + reward_zone_time + bottom_right_time);
    %all_percent_platform_zone = platform_zone_time/(top_left_time + platform_zone_time + reward_zone_time + bottom_right_time);
    all_percent_reward_zone = reward_zone_time/(platform_zone_time + reward_zone_time);
    all_percent_platform_zone = platform_zone_time/(platform_zone_time + reward_zone_time);

    %sound time
    sound_rows = rows;
    sound_columns = columns;
    sound_top_left = sound_time_spent(1:sound_rows, 1:sound_columns);
    sound_bottom_left = sound_time_spent((sound_rows + 1):end, 1:sound_columns);
    sound_top_right = sound_time_spent(1:sound_rows, (sound_columns + 1):end);
    sound_bottom_right = sound_time_spent((sound_rows + 1):end, (sound_columns + 1):end);

    if ~flip_arena
        sound_top_left_time = sum(sum(sound_top_left));
        sound_platform_zone_time = sum(sum(sound_top_right));
        sound_reward_zone_time = sum(sum(sound_bottom_left));
        sound_bottom_right_time = sum(sum(sound_bottom_right));
    else
        sound_top_left_time = sum(sum(sound_top_left));
        sound_reward_zone_time = sum(sum(sound_top_right));
        sound_platform_zone_time = sum(sum(sound_bottom_left));
        sound_bottom_right_time = sum(sum(sound_bottom_right));
        %CHANGED
        %sound_top_left_time = sum(sum(sound_bottom_right));
        %sound_platform_zone_time = sum(sum(sound_bottom_left));
        %sound_reward_zone_time = sum(sum(sound_top_right));
        %sound_bottom_right_time = sum(sum(sound_top_left));
    end
    %sound_percent_reward_zone = sound_reward_zone_time/(sound_top_left_time + sound_platform_zone_time + sound_reward_zone_time + sound_bottom_right_time);
    %sound_percent_platform_zone = sound_platform_zone_time/(sound_top_left_time + sound_platform_zone_time + sound_reward_zone_time + sound_bottom_right_time);
    sound_percent_reward_zone = sound_reward_zone_time/(sound_platform_zone_time + sound_reward_zone_time);
    sound_percent_platform_zone = sound_platform_zone_time/(sound_platform_zone_time + sound_reward_zone_time);
    
    %sound off time
    sound_off_top_left = sound_off_time_spent(1:sound_rows, 1:sound_columns);
    sound_off_bottom_left = sound_off_time_spent((sound_rows + 1):end, 1:sound_columns);
    sound_off_top_right = sound_off_time_spent(1:sound_rows, (sound_columns + 1):end);
    sound_off_bottom_right = sound_off_time_spent((sound_rows + 1):end, (sound_columns + 1):end);
    
    if ~flip_arena
        sound_off_top_left_time = sum(sum(sound_off_top_left));
        sound_off_platform_zone_time = sum(sum(sound_off_top_right));
        sound_off_reward_zone_time = sum(sum(sound_off_bottom_left));
        sound_off_bottom_right_time = sum(sum(sound_off_bottom_right));
    else
        sound_off_top_left_time = sum(sum(sound_off_top_left));
        sound_off_reward_zone_time = sum(sum(sound_off_top_right));
        sound_off_platform_zone_time = sum(sum(sound_off_bottom_left));
        sound_off_bottom_right_time = sum(sum(sound_off_bottom_right));
        %CHANGED
        %sound_off_top_left_time = sum(sum(sound_off_bottom_right));
        %sound_off_platform_zone_time = sum(sum(sound_off_bottom_left));
        %sound_off_reward_zone_time = sum(sum(sound_off_top_right));
        %sound_off_bottom_right_time = sum(sum(sound_off_top_left));
    end
    %sound_off_percent_reward_zone = sound_off_reward_zone_time/(sound_off_top_left_time + sound_off_platform_zone_time + sound_off_reward_zone_time + sound_off_bottom_right_time);
    %sound_off_percent_platform_zone = sound_off_platform_zone_time/(sound_off_top_left_time + sound_off_platform_zone_time + sound_off_reward_zone_time + sound_off_bottom_right_time);

    sound_off_percent_reward_zone = sound_off_reward_zone_time/(sound_off_platform_zone_time + sound_off_reward_zone_time + sound_off_bottom_right_time);
    sound_off_percent_platform_zone = sound_off_platform_zone_time/(sound_off_platform_zone_time + sound_off_reward_zone_time);
    
    
    save('Behavior/Sound on data', 'sound_reward_zone_time', 'sound_platform_zone_time', 'sound_top_left_time', 'sound_bottom_right_time', 'sound_percent_reward_zone', 'sound_percent_platform_zone');
    save('Behavior/Sound off data', 'sound_off_reward_zone_time', 'sound_off_platform_zone_time', 'sound_off_top_left_time', 'sound_off_bottom_right_time', 'sound_off_percent_reward_zone', 'sound_off_percent_platform_zone');
    save('Behavior/All data', 'reward_zone_time', 'platform_zone_time', 'top_left_time', 'bottom_right_time', 'all_percent_reward_zone', 'all_percent_platform_zone');

    %Graph those percents
    X = categorical({'All', 'Sound On', 'Sound Off'});
    Y = [all_percent_platform_zone, sound_percent_platform_zone, sound_off_percent_platform_zone];
    bar(X, Y);
    title('Percent time in platform zone');
    savefig('Behavior/platform_zone_percent.fig')

    %Graph those percents
    X = categorical({'All', 'Sound On', 'Sound Off'});
    Y = [all_percent_reward_zone, sound_percent_reward_zone, sound_off_percent_reward_zone];
    bar(X, Y);
    title('Percent time in reward zone');
    savefig('Behavior/reward_zone_percent.fig')

    new_time_filter = bandpass_disk_filter(new_time_spent, 20, 1000);
    sound_time_filter = bandpass_disk_filter(sound_time_spent, 20, 1000);
    sound_off_time_filter = bandpass_disk_filter(sound_off_time_spent, 20, 1000);

    new_time_spent_sum = sum(sum(new_time_spent));
    sound_time_spent_sum = sum(sum(sound_time_spent));
    sound_off_time_spent_sum = sum(sum(sound_off_time_spent));

    sound_correction = new_time_spent_sum/sound_time_spent_sum;
    sound_off_correction = new_time_spent_sum/sound_off_time_spent_sum;

    corrected_sound_time_filter = bandpass_disk_filter(sound_time_spent .* sound_correction, 20, 1000);
    corrected_sound_off_time_filter = bandpass_disk_filter(sound_off_time_spent .* sound_off_correction, 20, 1000);

    %% plot and save data
    %%%%% all graphs
    h1 = figure; 
    set(gcf,'units','points','position',[400,10,900,400])
    subplot(2,3,1); 
    %full vertical, half horizontal 1

    %new_J = imresize(J, [400 300]); 
    new_J = J;
    imagesc(new_J,[0 mean2(new_J)+3*std2(new_J) - 100]); 
    axis('square')
    %image(J);
    colormap(gca,'gray')
    hold on; 
    %plot(data(:,1),data(:,2),'r');
    plot(smoothed_track_data(:,1),smoothed_track_data(:,2),'c');
    %plot(smoothed_track_data(:,1),smoothed_track_data(:,2),'r');
    box off; axis off; 
    title('all trajectory');

    subplot(2,3,4);
    %%% TODO fix the heatmap
    %imagesc(new_time_filter(1:400, 100:400),[0 .04]); 
    temp;
    imagesc(new_time_filter,[0 .4]);
    % colormap(gca,'jet');
    axis('square')
    box off; axis off; 
    title('all heatmap');

    %%%%%%%%%%% sound graphs
    subplot(2,3,2); 
    %full vertical, half horizontal 1
    imagesc(J,[0 mean2(J)+3*std2(J)]); 
    colormap(gca,'gray')
    
    hold on; 
    %plot(data(:,1),data(:,2),'r');
    plot(sound_on_data(:,1),sound_on_data(:,2),'c');
    axis('square')
    box off; axis off; 
    title('sound on trajectory');

    subplot(2,3,5);
    %%% TODO fix the heatmap
    imagesc(corrected_sound_time_filter,[0 .4]); 
    axis('square')
    % colormap(gca,'jet');
    box off; axis off; 
    title('sound on heatmap');

    %%%%%%%%%%%%%
    %Sound off graphs
    subplot(2,3,3); 
    %full vertical, half horizontal 1
    imagesc(J,[0 mean2(J)+3*std2(J)]); 
    colormap(gca,'gray')
    hold on; 
    %plot(data(:,1),data(:,2),'r');
    plot(sound_off_data(:,1),sound_off_data(:,2),'c');
    axis('square');
    box off; axis off; 
    title('sound off trajectory');

    subplot(2,3,6);
    %%% TODO fix the heatmap
    imagesc(corrected_sound_off_time_filter,[0 .5]); 
    axis('square')
    % colormap(gca,'jet');
    box off; axis off; 
    title('sound off heatmap');

    %% TODO normalize by time in each condition. Sound is too weak DONE
    %% TODO output graph with percent time spent in reward zone vs platform zone in all conditions
    %% save a background image so you don't have to keep so many video files
    save(file_nm, 'mouse_parameter');
    saveas(h1,file_nm)
    close all
    %% Find likelihood on platform during sound
    [platform_matrix, sound_data_time_adjusted] = find_sound_platform_position(sound_on_data, J);
    
    starts = find(sound_data_time_adjusted(:,3) == 0);
    
    segmented_platform_matrix = [];
    for p = 1:size(starts, 1)
        prev = starts(p, 1);
        segment = transpose(platform_matrix(prev:(prev + 599), 1));
        segmented_platform_matrix(p, 1: 600) = segment;
        
    end
    temp_axis = 0:0.0333:20;
    temp_axis = temp_axis(2:end);
    platform_during_sound = smooth(mean(segmented_platform_matrix));
    plot(temp_axis, platform_during_sound, 'LineWidth', 3);
    title('Likelihood of being on the platform during the sound', 'FontSize', 20)
    xlabel('Seconds from sound onset', 'FontSize', 20);
    ylim(gca, [0,1]);
    savefig('Behavior/platform_likelihood.fig')
    save('Behavior/platform_likelihood.mat', 'temp_axis', 'platform_during_sound');
    
    %% Find likelihood on platform during 20 s pre sound
    sound_on_times_temp = sound_on_times;
    times = track_data(:,3);
    pre_sound_inds = [];
    full_matrix_inds = [];
    trial_num = 1;
    for x = 1:size(times, 1)
        time = times(x);
        if isempty(sound_on_times_temp)
            break
        end
        if time > (sound_on_times_temp(1) - 40)
            if time < sound_on_times_temp(1)
                %If in 20 seconds pre sound
                pre_sound_inds = horzcat(pre_sound_inds, x);
            else
                %Once it hits sound on
                sound_on_times_temp(1) = [];
                num = size(pre_sound_inds, 2);
                if num < 1200
                    dif = 1200 - num;
                    to_add = ones(1, dif);
                    pre_sound_inds = horzcat(to_add, pre_sound_inds);
                end
                full_matrix_inds = vertcat(full_matrix_inds, pre_sound_inds(:, 1:1200));
                pre_sound_inds = [];
            end
        end
    end
    
    full_platform_matrix_sound_off = [];
    for g = 1:size(full_matrix_inds, 1)
        [platform_matrix_off, sound_off_data_time_adjusted] = find_sound_platform_position(track_data(full_matrix_inds(g,:), :), J);
        full_platform_matrix_sound_off = horzcat(full_platform_matrix_sound_off, platform_matrix_off(:, 1));
    end
    
    temp_axis = 0:0.033333:40;
    temp_axis = temp_axis(2:end);
    avg_full_platform_matrix_sound_off = smooth(mean(full_platform_matrix_sound_off, 2));
    
    %% Find likelihood on platform during 20 s post sound 
    sound_on_times_temp = sound_on_times;
    post_sound_inds = [];
    post_matrix_inds = [];
    for x = 1:size(times, 1)
        time = times(x);
        if isempty(sound_on_times_temp)
            break
        end
        if time > (sound_on_times_temp(1) + 20)
            if time < (sound_on_times_temp(1) + 60)
                %If in 20 seconds post sound
                post_sound_inds = horzcat(post_sound_inds, x);
            else
                %Once it hits sound on
                sound_on_times_temp(1) = [];
                num = size(post_sound_inds, 2);
                if num < 1200
                    last = post_sound_inds(end);
                    dif = 1200 - num;
                    to_add = (last + 1):(last + dif);
                    post_sound_inds = horzcat(post_sound_inds, to_add);
                end
                post_matrix_inds = vertcat(post_matrix_inds, post_sound_inds(:, 1:1200));
                post_sound_inds = [];
            end
        end
    end
    
    full_platform_matrix_post_sound = [];
    for g = 1:size(post_matrix_inds, 1)
        [platform_matrix_post, ~] = find_sound_platform_position(track_data(post_matrix_inds(g,:), :), J);
        full_platform_matrix_post_sound = horzcat(full_platform_matrix_post_sound, platform_matrix_post(:, 1));
    end
    
    temp_axis = 0:0.033333:100;
    temp_axis = temp_axis(2:end);
    avg_full_platform_matrix_post_sound = smooth(mean(full_platform_matrix_post_sound, 2));
    
    %% Plot full trial
    
    
    full_trial = vertcat(avg_full_platform_matrix_sound_off, platform_during_sound, avg_full_platform_matrix_post_sound);
    
    plot(temp_axis, full_trial, 'LineWidth', 3);
    title('Full trial', 'FontSize', 20)
    ylabel('P(Platform)', 'FontSize', 20);
    xlabel('Seconds', 'FontSize', 20);
    hold on 
    s = plot([40, 40], [0, 1], 'r');
    plot([60, 60], [0, 1], 'r');
    legend(s, 'sound');
    ylim(gca, [0,1]);
    savefig('Behavior/platform_likelihood_full.fig')
    save('Behavior/platform_likelihood_full.mat', 'temp_axis', 'full_trial');
    
end



if record_tracking_video
    close all
    reader = VideoReader(horzcat(file_nm, file_type));
    writer = VideoWriter('video_with_arrow', 'MPEG-4');
    writer.FrameRate = reader.FrameRate;
    start_of_bpod_frame = start_of_bpod_sec * reader.FrameRate;
    open(writer);
    p = 1;
    max_size = size(smoothed_track_data, 1);
    if record_sample
        max_size = 7200;
        FrameNumber = 7200;
    end
    hw = waitbar(0,'Animal tracking ...');
    w = 1;
    while and(hasFrame(reader), p <= max_size)
        tic; % start of clock
        img = readFrame(reader);
        
        if p < start_of_bpod_frame + 5
            p = p + 1;
            continue
        end
        I2 = im2double(img);
        I2 = imcrop(I2, rect2);
        %centroid_y = round(smoothed_track_data(w, 1));
        %centroid_x = round(smoothed_track_data(w, 2));
        centroid_y = round(track_data(w, 1));
        centroid_x = round(track_data(w, 2));
        w = w + 1;
        if centroid_x < 4
            centroid_x = 4;
        end
        if centroid_y < 4
            centroid_y = 4;
        end
        %nose_x = round(nose_points(p, 1));
        %nose_y = round(nose_points(p, 2));
        %image8Bit = uint8(255 * mat2gray(floatingPointImage));
        I2((centroid_x - 3):(centroid_x + 3), (centroid_y - 3):(centroid_y + 3), 1) = 1;
        I2((centroid_x - 3):(centroid_x + 3), (centroid_y - 3):(centroid_y + 3), 2) = 1;
        I2((centroid_x - 3):(centroid_x + 3), (centroid_y - 3):(centroid_y + 3), 3) = 1;

        %plot orientation
        %I2((centroid_x + nose_x - 3):(centroid_x + nose_x + 3), (centroid_y + nose_y - 3):(centroid_y + nose_y + 3), 1) = 0;
        %I2((centroid_x + nose_x - 3):(centroid_x + nose_x + 3), (centroid_y + nose_y - 3):(centroid_y + nose_y + 3), 2) = 1;
        %I2((centroid_x + nose_x - 3):(centroid_x + nose_x + 3), (centroid_y + nose_y - 3):(centroid_y + nose_y + 3), 3) = 1;
        imshow(I2);
        writeVideo(writer,I2);

        if mod(p, 100) == 0
            e_time = toc; % time per loop
            e_time_all = e_time*FrameNumber;
            remaining_time = e_time_all - e_time*p;
            waitbar(p/FrameNumber,hw,['remaining time ',num2str(remaining_time,'%2.0f'),'s']);
        end

        p = p + 1;
    end
    close(hw)
    close(writer);
end
end
function [photometry_start_time, start_times] = get_start_times(file, analog_key_word)
    file = horzcat(file, '/'); 
    start_file = dir(horzcat(file, '*', analog_key_word, '*.csv'));
    if ~isempty(start_file)
        photometry_start = dlmread(horzcat(file, start_file.name));
    else
        photometry_start = [0];
    end
    photometry_start_time = photometry_start(1);
    if 1
        %Convert analog in file to seconds
        photometry_start_clone = photometry_start;
        photometry_start_clone = (photometry_start_clone - photometry_start(1,1))/1000;
        if photometry_start_clone == 0
            error('Change the analog_key_word to something in the name of your file.');
        end
        if photometry_start_clone(2) - photometry_start_clone(1) > 1
            photometry_start_clone(1) = [];
            photometry_start(1,:) = [];
            photometry_start_clone = photometry_start_clone - photometry_start_clone(1);
        end
        photometry_start_clone = photometry_start_clone';

        prev = 1;
        cur = 1;
        alt_start_times = [photometry_start_clone(1)];
    else
        %Convert centroid and analog file to seconds
        photometry_start_clone = photometry_start;
        photometry_start_clone = (photometry_start_clone - photometry_start(1,1))/1000;
        if photometry_start_clone(2) - photometry_start_clone(1) > 1
            photometry_start_clone(1) = [];
            photometry_start(1,:) = [];
            photometry_start_clone = photometry_start_clone - photometry_start_clone(1);
        end
        photometry_start_clone = photometry_start_clone';

        prev = 1;
        cur = 1;
        alt_start_times = [photometry_start_clone(1)];
    end

    for s = 2:size(photometry_start_clone,2)
        dif = photometry_start_clone(s) - photometry_start_clone(prev);
        if dif > 5
            alt_start_times = horzcat(alt_start_times, photometry_start_clone(s));
        end
        prev = prev + 1;
    end

    %% Use this to keep track of matlab events
    %set 0 point for imaging beginning. This coincides with trial 1.
    startRef = photometry_start(1);

    %convert time to relative time (s)
    start_times = alt_start_times;
end

function [platform_matrix, sound_data_time_adjusted] = find_sound_platform_position(tracking_data, background_image)
%input- tracking data N x 3 (x,y,time from onset of sound) output- matrix N
%x 2 (boolean, timestamp)

    sound_data_time_adjusted = tracking_data;
    prev = 0;
    index = 1;
    for y = 1:size(tracking_data, 1)
        cur = tracking_data(y, 3);
        if prev == 0
            prev = cur;
            continue
        end
        if cur > prev + 20
            sound_data_time_adjusted(index:(y-1), 3) = sound_data_time_adjusted(index:(y-1), 3) - sound_data_time_adjusted(index, 3);
            index = y;
            prev = cur;
        end
    end

    platform_matrix = zeros(size(tracking_data, 1), 2);
    for i = 1:size(tracking_data, 1)
        vert_num = round(tracking_data(i, 2));
        horz_num = round(tracking_data(i, 1));
        time = tracking_data(i, 3);
        
        
        %I2 = im2double(background_image);
        %I2((abs(vert_num - 3)):(vert_num + 3), (abs(horz_num - 3)):(horz_num + 3), 1) = 0;
        %imshow(I2);
        
        % Find if the mouse is on the platform
        if vert_num > size(background_image, 1)/2 + 10
            %if in bottom half
            if horz_num < (size(background_image, 2)/2 - 10)
                %if in left half
                %disp('In platform zone');
                platform_matrix(i, 1:2) = [1, time];
            else
                %disp('Not in platform zone');
                platform_matrix(i, 1:2) = [0, time];
            end
        else
            %disp('Not in platform zone');
            platform_matrix(i, 1:2) = [0, time];
        end
    end
end
