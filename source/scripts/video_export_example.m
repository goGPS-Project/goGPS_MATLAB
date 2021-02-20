num_frame = 120; % number of frames of the animation
f_handle = gcf; % get the handle of the figure

% Set up video export ---------------------------------------------------------

fprintf('Exporting frame %5d/%5d', 0, num_frame);
if ismac() || ispc()
    % Better compression on Mac > 10.7 and Win > 7
    video_out = VideoWriter(fullfile(Core.getState.getOutDir, ['AniMap' par_str_short '.mp4']), 'MPEG-4');
else
    % Linux doesn't have mp4 compression avaiable
    video_out = VideoWriter(fullfile(Core.getState.getOutDir, ['AniMap' par_str_short '.avi']));
end
video_out.FrameRate = 30;
video_out.Quality = 91;
open(video_out);

% Start animation -------------------------------------------------------------
for i = 1 : num_frame
    fprintf('%s%5d/%5d', char(8 * ones(1,11)), i, num_frame);
    
    % -----------------------
    %  UPDATE YOUR PLOT HERE
    % -----------------------
    
    frame = getframe(f_handle); % get the frime of figure f
    writeVideo(video_out, frame);
end
close(video_out);
