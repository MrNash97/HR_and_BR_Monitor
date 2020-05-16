%cam = videoinput('winvideo', 2, 'RGB24_640x480');
cam = videoinput('winvideo', 1, 'MJPG_640x480');
%cam = videoinput('winvideo', 2, 'MJPG_640x480');
src = getselectedsource(cam);
src.FrameRate = '30.0000';
%src.ExposureMode = 'manual';
%src.WhiteBalanceMode = 'manual';
cam.FramesPerTrigger = inf;
% Skip the first few frames the device provides
% before logging data.
vidobj.TriggerFrameDelay = 5;

videoFrame = getsnapshot(cam);
frameSize = size(videoFrame);

videoPlayer = vision.VideoPlayer('Position', [100 0 [frameSize(2), frameSize(1)]+30]);

runLoop = true;
numPts = 0;
frameCount = 0;
start(cam);
tic;
while(cam.FramesAcquired<500)
    videoFrame=getdata(cam,1,'uint8');
    step(videoPlayer, videoFrame);
end
t=toc;
stop(cam);
release(videoPlayer);
flushdata(cam);