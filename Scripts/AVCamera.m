% create a video obj to get data
vidobj = videoinput('gentl', 1, 'Mono8'); 
vidobj.TriggerRepeat=inf;
vidobj.FramesPerTrigger=1;
% set trigger mode to manual, then we can start and stop the
% camera by our self
triggerconfig(vidobj, 'manual');

% Get the video source
src = getselectedsource(vidobj); 


% 
% % set the Binning of the camera
% src.BinningHorizontal = 2;
% src.BinningVertical = 2;

% start video aquisition
start(vidobj);
% get a pic
trigger(vidobj);
%snapshot = getsnapshot(vidobj); 
snapshot=getdata(vidobj);

% adjust to the size of our pic
im = image(zeros(size(snapshot),'uint8'));

% set Preview as image mode
axis('image');


for i = 1:3
    % update Exposure Time and Gain

    src.ExposureTime = 112.181;
    src.Gain = 1;

    % grab a new figure from camera
    snapshot = getsnapshot(vidobj);


    % update Preview figure
    imagesc(snapshot,[0 255]);

    % Pause to get stop
    pause(0.1);

end

% we must stop this video object by ourself, otherwise problem
% will come
stop(vidobj);

% delete the video object, so other thread can use this camera
delete(vidobj);
% clear;