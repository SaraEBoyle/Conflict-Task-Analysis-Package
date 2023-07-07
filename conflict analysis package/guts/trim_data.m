function [data, auto_flo] = trim_data(data, fiber, subtract_autofluorescence)
    %% Cut data from beginning and/or end and subtract baseline fluorescence
    %Plot pre-bleach corrected data.
    if subtract_autofluorescence
        plot(data(:,1), data(:,fiber + 1));
        figure;
        if size(data, 1) > 8000
            plot(data(1:8000,1), data(1:8000,fiber + 1));
        else
            plot(data(1:size(data, 1),1), data(1:size(data, 1),fiber + 1));
        end
        title('First 800 seconds. Choose autofluorescence window. Click at beginning, then at end')
        [x_in,y_in]=ginput(2);
        disp(horzcat('Window is from ', num2str(x_in(1)), ' to ', num2str(x_in(2))))
        %disp(y_in)
        close all
        is_good = input('Did you choose? Press enter to continue, e to exit: ', 's');
        if is_good == 'e'
            disp('User exited. Closed program.')
            close all
            return
        end
    end
    %There's weird frame dropping
    if subtract_autofluorescence == 1
        auto = data(find(data(:,1) > x_in(1)), :);
        auto = auto(find(auto(:,1) < x_in(2)), fiber + 1);
        auto_flo = mean(auto);
        if isnan(auto)
            error('You chose autofluorescence times wrong.')
        end
        %auto_flo = 5.28;
    else 
        auto_flo = 0;
    end

    %Cut out beginning of file
    data = data(find(data(:,1) > -12), :);

    %Check if this looks good
    plot(data(:,1), data(:,fiber + 1));
    is_good = input('Does this look right? Press enter to continue without cutting any data, c to cut, e to exit: ', 's');
    if is_good == 'e'
        disp('User exited. Closed program.')
        close all
        return
    elseif is_good == 'c'
        disp('Select interval to analyze. ');
        [x_in,y_in] = ginput(2);

        %Cut anything funky from the end
        data = data(find(data(:,1) > x_in(1)), :);
        data = data(find(data(:,1) < x_in(2)), :);
    end
    close all;
end