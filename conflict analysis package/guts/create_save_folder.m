function [] = create_save_folder(fiber_name, alignment, file)
    %Sets up folders to save analyzed data, cds to new folder. Reads the bpod file in.

    %% Set up save folders
    %Create the fiber specific folders
    if ~exist(horzcat(file,fiber_name), 'dir')
        mkdir(fiber_name);
        cd(fiber_name);
    elseif exist(horzcat(file,fiber_name), 'dir')
        disp(horzcat(fiber_name, ' trials directory already exists.'))
        current_dir = pwd;
        if ~strcmp(horzcat(current_dir, '/'), horzcat(file, fiber_name, '/'))
            cd(fiber_name);
        else
            disp('You are in the correct directory');
        end
    end
    

    %Create the folder that specifies trial type 
    alignment_folder = alignment;
    if ~exist(alignment_folder, 'dir')
            mkdir(alignment_folder);
            cd(alignment_folder);
        elseif exist(alignment_folder, 'dir')
            disp(horzcat(alignment_folder, ' trials directory already exists.'))
            current_dir = pwd;
            if ~strcmp(horzcat(current_dir, '/'), horzcat(file, fiber_name, '/', alignment_folder, '/'))
                cd(alignment_folder);
            else
                disp('You are in the correct directory');
            end
    end
    close all
end