function [gcamp, iso] = de_interleave(trigger1, trigger2, file, fiber_name, data, fiber, auto_flo, subtract_autofluorescence, signal_v)
%% de-interleave and subtract auto fluorescence
    %Use this if you alternate channels. Can handle 2 channels
    if ~and(trigger2, and(exist(horzcat(file, fiber_name, '/', 'signal.mat')) == 2, exist(horzcat(file, fiber_name, '/', 'control.mat')) == 2))
        if trigger1 == 1
            gcamp = data(:, [1 fiber + 1]);
            avg_F = mean(gcamp(:,2));
            temp_fit = fit(gcamp(:,1),gcamp(:,2),'exp2'); %find bi-exponential best fit
            close all
            %fitted is the line that forms the threshold to split channels
            adjust_threshold = 0;
            fitted = temp_fit(gcamp(:,1)) + adjust_threshold;
            for tries = 1:5
                figure
                plot(gcamp(:,1),gcamp(:,2)) %plot data and fit on same graph to make sure fit looks ok
                hold on
                fitted = fitted + adjust_threshold;

                plot(gcamp(:,1), fitted);
                title('Threshold to divide channels')
                is_good = input('Check whether the channel separation is good. Press enter if yes, n to adjust threshold, e to exit. ', 's');

                if strcmp(is_good, 'e')
                    close all
                    disp('User exited.')
                    return
                elseif strcmp(is_good, 'n')
                    new_thresh = input('Enter value to add to threshold line to adjust: (+ or -). e to exit ', 's');
                    if strcmp(new_thresh, 'e')
                        break
                    else
                        adjust_threshold = str2num(new_thresh);
                        hold off 
                        close all
                    end
                    continue
                else
                    break
                end
            end
            if tries == 5
                error('Are you sure your threshold is ok?');
            end
            close all;
            %plot(temp_fit);
            original = gcamp(:,2);
            %fitted = temp_fit(gcamp(:,1));
            hold off
            close all
            data_1 = [];
            data_2 = [];
            for x = 1:size(gcamp, 1)
                time = gcamp(x, 1);
                sig = gcamp(x, 2);
                %set fit as boundary threshold to separate channels
                thresh = fitted(x);
                if gcamp(x, 2) > thresh
                    data_1 = vertcat(data_1, [time, sig]);         
                else
                    data_2 = vertcat(data_2, [time, sig]);           
                end
            end
            figure;
            subplot(2,1,1)
            plot(data_1(:,1), data_1(:,2));
            title('Channel 1')
            subplot(2,1,2)
            plot(data_2(:,1), data_2(:,2));
            title('Channel 2')
            if signal_v == 1
                which_channel = input('Is the signal from channel 1 or 2? e to exit, answer either 1 or 2: ', 's');
            else
                which_channel = input('Is the control from channel 1 or 2? e to exit, answer either 1 or 2: ', 's');
            end
            

            if strcmp(which_channel, 'e')
                close all
                return
            end
            if trigger2 == 1
                channel_type = input('Control channel (c) or signal (s)? ','s');
            end
            close all;
            if str2num(which_channel) == 1
                if trigger2 == 0
                    gcamp = data_1;
                    iso = data_2;
                else
                    if strcmp(channel_type, 's')
                        gcamp = data_1;
                        iso = data_2;
                    elseif strcmp(channel_type, 'c')
                        gcamp = data_2;
                        iso = data_1;
                    end

                end
            elseif str2num(which_channel) == 2
                if trigger2 == 0
                    gcamp = data_2;
                    iso = data_1;
                else
                    if strcmp(channel_type, 's')
                        gcamp = data_2;
                        iso = data_1;
                    elseif strcmp(channel_type, 'c')
                        gcamp = data_1;
                        iso = data_2;
                    end
                end

            end
        end
    end

    if trigger2 == 1
        %Specify which channel you're looking at and save
        if and(exist(horzcat(file, fiber_name, '/', 'signal.mat')) == 2, exist(horzcat(file, fiber_name, '/', 'control.mat')) == 2)
            load(horzcat(file, fiber_name, '/', 'control.mat'))
            load(horzcat(file, fiber_name, '/', 'signal.mat'))

        elseif strcmp(channel_type, 's')
            auto_flo_s = auto_flo;
            if ~exist(horzcat(file, fiber_name, '/', 'signal.mat'), 'var')   
                save(horzcat(file, fiber_name, '/', 'signal.mat'), 'gcamp', 'auto_flo_s')
            end
            disp('Add a control channel. Exited.');
            cd(file)
            return;

        elseif strcmp(channel_type, 'c')
            auto_flo_c = auto_flo;
            if ~exist(horzcat(file, fiber_name, '/', 'control.mat'), 'var')
                save(horzcat(file, fiber_name, '/', 'control.mat'), 'iso', 'auto_flo_c')
            end
            disp('Add a signal channel if you have not yet. Exited.');
            cd(file)
            return;
        end

        %If you have saved control and signal, load them as iso_trials and
        %trials
    elseif trigger1 == 1
        auto_flo_s = auto_flo;
        auto_flo_c = auto_flo;
        %gcamp = data(:, [1 fiber + 1]);
    else
        auto_flo_s = auto_flo;
        auto_flo_c = auto_flo;
        gcamp = data(:, [1 fiber + 1]);
    end
    
    %% Subtract autofluorescence
    if subtract_autofluorescence == 1
        gcamp(:,2) = gcamp(:,2) - auto_flo_s;
        if or(trigger1, trigger2)
            iso(:,2) = iso(:,2) - auto_flo_c;
        end
    else
        auto_flo = 0;
    end
end