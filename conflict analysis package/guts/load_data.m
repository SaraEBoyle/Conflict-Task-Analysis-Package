function[start_times, data, SessionData, new_data] = load_data(bpod_key_word, signal_key_word, analog_key_word, dir_nm, photometry_data)
    %% select EthoVision data from excel
    % clear nam
    %% Load in data
    %Retrieve the bpod behavior file
    file = dir_nm;
    tem = dir(horzcat(file, '*', bpod_key_word, '*.mat')); 
    if isempty(tem)
        tem = dir(horzcat(file, '*', bpod_key_word, '*.mat'));
        if isempty(tem)
            tem = dir(horzcat(file, '*behavior*.mat'));
        end
    end
    if isempty(tem)
        error('set bpod_key_word to a phrase that appears in the bpod file')
    end

    name = tem(1).name;
    folder = tem(1).folder;
    full_path = horzcat(folder, '/', name);
    disp(name);
    load(full_path);
    
    if photometry_data
    %Load photometry data
        signal_file = dir(horzcat(file, '*', signal_key_word, '*.csv'));
        data = dlmread(horzcat(file, signal_file.name));
    end
    
    centroid_file = dir(horzcat(file, '*', 'centroid', '*.csv'));
    new_data = dlmread(horzcat(file, centroid_file.name));
    
    start_file = dir(horzcat(file, '*', analog_key_word, '*.csv'));
    if ~isempty(start_file)
        photometry_start = dlmread(horzcat(file, start_file.name));
    else
        photometry_start = [0];
    end
    
    if photometry_data
        %Convert analog in file to seconds
        photometry_start_clone = photometry_start;
        photometry_start_clone = (photometry_start_clone - photometry_start(1,1))/1000;
        data(:,1) = (data(:,1) - photometry_start(1,1))/1000;
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
        data(:,1) = (data(:,1) - photometry_start(1,1))/1000;
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