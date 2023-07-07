
%for trial 16, 33.2167 minutes is sound, or 1993 secs. Do 20 seconds pre
%and 40 post, or 1973-2033, or 32:43-33:

%Sound - 33:13, start: 32:53, end: 33:53
%gcamp = raw_signal_dFFs{1};
%iso = raw_control_dFFs{1};
%%%%%%%%%%%% winner winner. 3rd bought
sig_y = gcamp(:, 2).*100;
cont_y = iso(:, 2).*100;
signal = horzcat(gcamp(:, 1), sig_y);
control = horzcat(iso(:, 1), cont_y);
fps = 10;

%% 5 min 8.0580, start is 4 min 48.0600

%For black, the lims are -1 and 2
%For purple 176.007846400000 before bpod start
% [377.115545600005] trial 3. Corrected for movie it is  437.0590

movie_offset = 96.5815;
event_start = 1896.455462400004; %in seconds
time_before = 20;
time_after = 40;
movie_time = event_start + movie_offset;
disp(horzcat('Sound time: ', num2str(floor(movie_time/60)), 'min, ', num2str((movie_time/60 - floor(movie_time/60)) * 60), ' sec'));
disp(horzcat('Start time: ', num2str(floor((movie_time - time_before)/60)), 'min, ', num2str(((movie_time - time_before)/60 - floor((movie_time - time_before)/60)) * 60), ' sec'));
disp(horzcat('End time: ', num2str(floor((movie_time + time_after)/60)), 'min, ', num2str(((movie_time + time_after)/60 - floor((movie_time + time_after)/60)) * 60), ' sec'));
movie_start = event_start - time_before;
movie_end = event_start + time_after;
%3071 is 325,  3181 frame is 336 seconds, third wc bought
% drop is at 401, go until 415. 3971
biggers = find(signal(:, 1) >= movie_start);
start_frame = biggers(1);
bigger = find(signal(:, 1) >= movie_end);
end_frame = bigger(1);
signal(:, 1) = signal(:, 1) - signal(start_frame, 1);
control(:, 1) = control(:, 1) - control(start_frame, 1);

frames_per_movie = end_frame - start_frame;
%% 0 s to [107.782323199995]
for t = 1:frames_per_movie
   %fplot(@(x) sin(x*50/t),[0,2*pi]);  % plot
   plot(signal((start_frame):(start_frame + t), 1), signal((start_frame):(start_frame + t), 2), 'color', 'g', 'LineWidth', 2);
   hold on
   xlabel('Seconds');
   ylabel('dF/F (%)')
   plot(control((start_frame):(start_frame + t), 1), control((start_frame):(start_frame + t), 2), 'color', 'm', 'LineWidth', 2);
   ylim([-1,2]);                      % guarantee consistent height
   %xlim([signal(1 + t, 1), signal(1000 + t, 1)]);
   xlim([signal(start_frame, 1), signal(start_frame + frames_per_movie, 1)]);
   F(t) = getframe;                   % capture it
   hold off
end

%% fourth bought
%4968 frame 515 s
%500 s 4819
% 4669 4.850538495999947e+02
% [488.157977600001] s [578.270169599995] s
%8 7.8, 9 38.27

writerObj = VideoWriter('trial_16.avi');
writerObj.FrameRate = 10;
open(writerObj);
writeVideo(writerObj, F)
close(writerObj);
