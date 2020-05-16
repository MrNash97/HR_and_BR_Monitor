%% Face Detection and Tracking Using Live Video Acquisition
% This example shows how to automatically detect and track a face in a live
% video stream, using the KLT algorithm.   
% 
%   Copyright 2014 The MathWorks, Inc.

%% Overview
% Object detection and tracking are important in many computer vision
% applications including activity recognition, automotive safety, and
% surveillance.  In this example you will develop a simple system for
% tracking a single face in a live video stream captured by a webcam.
% MATLAB provides webcam support through a Hardware Support Package,
% which you will need to download and install in order to run this example. 
% The support package is available via the 
% <matlab:supportPackageInstaller Support Package Installer>.
%
% The face tracking system in this example can be in one of two modes:
% detection or tracking. In the detection mode you can use a
% |vision.CascadeObjectDetector| object to detect a face in the current
% frame. If a face is detected, then you must detect corner points on the 
% face, initialize a |vision.PointTracker| object, and then switch to the 
% tracking mode. 
%
% In the tracking mode, you must track the points using the point tracker.
% As you track the points, some of them will be lost because of occlusion. 
% If the number of points being tracked falls below a threshold, that means
% that the face is no longer being tracked. You must then switch back to the
% detection mode to try to re-acquire the face.

%% Setup
% Create objects for detecting faces, tracking points, acquiring and
% displaying video frames.

% Create the face detector object.
%faceDetector = vision.CascadeObjectDetector();
faceDetector = vision.CascadeObjectDetector('FrontalFaceLBP');
%faceDetector = vision.CascadeObjectDetector('ClassificationModel','RightEye');

% Create the point tracker object.
pointTracker = vision.PointTracker('MaxBidirectionalError', 1);
%pointTracker = vision.PointTracker('NumPyramidLevels', 4);

% Create the webcam object.
%% build-in webcam
cam = webcam('EasyCamera');
cam.Resolution = '640x480';
%cam.Resolution = '320x240';
cam.WhiteBalanceMode = 'auto';
%cam.WhiteBalance = 4500;
%% External webcam
%% Logitech
% cam = webcam('Logitech');
% cam.Resolution = '640x480';
% cam.ExposureMode = 'auto';
% cam.Exposure = -4;
%%
% Capture one frame to get its size.

%videoFrame = step(cam);
videoFrame = snapshot(cam);
frameSize = size(videoFrame);
[height,width,~] = size(videoFrame);
% Create the video player object. 
videoPlayer = vision.VideoPlayer('Position', [100 0 [frameSize(2), frameSize(1)]+30]);

%% Detection and Tracking
% Capture and process video frames from the webcam in a loop to detect and
% track a face. The loop will run for 400 frames or until the video player
% window is closed.
trig = 0;
runLoop = true;
numPts = 0;
frameCount = 0;
t_tot = 0;
totalFrame = 500;
t_hist=zeros(1,totalFrame);
t_tot_hist=zeros(1,totalFrame);
rawColorSignal=zeros(3,totalFrame);
%% main program
while runLoop && frameCount < totalFrame
    if trig < 10
       videoFrame = snapshot(cam);
       videoFrame = im2double(videoFrame); %convert from 0-255 (uint8) to 0-1 (double) for masking purpose
       videoFrameGray = rgb2gray(videoFrame);
       trig = trig + 1;
    if numPts < 15
        bbox = faceDetector.step(videoFrameGray);
        if ~isempty(bbox)
            points = detectMinEigenFeatures(videoFrameGray, 'ROI', bbox(1, :),'MinQuality',0.01);
            xyPoints = points.Location;
            numPts = size(xyPoints,1);
            release(pointTracker);
            initialize(pointTracker, xyPoints, videoFrameGray);
            oldPoints = xyPoints;
            bboxPoints = bbox2points(bbox(1, :));  
            bboxPolygon = reshape(bboxPoints', 1, []);
        end
    else
        [xyPoints, isFound] = step(pointTracker, videoFrameGray);
        visiblePoints = xyPoints(isFound, :);
        oldInliers = oldPoints(isFound, :);
                
        numPts = size(visiblePoints, 1);       
        
        if numPts >= 15
            [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
                oldInliers, visiblePoints, 'similarity', 'MaxDistance', 3);            
            bboxPoints = transformPointsForward(xform, bboxPoints);
            bboxPolygon = reshape(bboxPoints', 1, []);            
            oldPoints = visiblePoints;
            setPoints(pointTracker, oldPoints);
        end
    end
    end
    tic;
    % Get the next frame.
    %videoFrame = step(cam);
    videoFrame = snapshot(cam);
    videoFrame = im2double(videoFrame); %convert from 0-255 (uint8) to 0-1 (double) for masking purpose
    videoFrameGray = rgb2gray(videoFrame);
    frameCount = frameCount + 1;

    if numPts < 15
        % Detection mode.
        bbox = faceDetector.step(videoFrameGray);
        
        if ~isempty(bbox)
            % Find corner points inside the detected region.
            
            points = detectMinEigenFeatures(videoFrameGray, 'ROI', bbox(1, :),'MinQuality',0.01);
            %points = detectHarrisFeatures(videoFrameGray, 'ROI', bbox(1, :));
            
            % Re-initialize the point tracker.
            xyPoints = points.Location;
            numPts = size(xyPoints,1);
            release(pointTracker);
            initialize(pointTracker, xyPoints, videoFrameGray);
            
            % Save a copy of the points.
            oldPoints = xyPoints;
            
            % Convert the rectangle represented as [x, y, w, h] into an
            % M-by-2 matrix of [x,y] coordinates of the four corners. This
            % is needed to be able to transform the bounding box to display
            % the orientation of the face.
            bboxPoints = bbox2points(bbox(1, :));  
            
            % Convert the box corners into the [x1 y1 x2 y2 x3 y3 x4 y4] 
            % format required by insertShape.
            bboxPolygon = reshape(bboxPoints', 1, []);
            
            % Display a bounding box around the detected face.
            videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, 'LineWidth', 3);
            
            % Display detected corners.
            videoFrame = insertMarker(videoFrame, xyPoints, '+', 'Color', 'white');
        end
        
    else
        % Tracking mode.
        [xyPoints, isFound] = step(pointTracker, videoFrameGray);
        visiblePoints = xyPoints(isFound, :);
        oldInliers = oldPoints(isFound, :);       
        numPts = size(visiblePoints, 1);       
        
        if numPts >= 15
            % Estimate the geometric transformation between the old points
            % and the new points.
            [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
                oldInliers, visiblePoints, 'similarity', 'MaxDistance', 3);            
            
            % Apply the transformation to the bounding box.
            bboxPoints = transformPointsForward(xform, bboxPoints);
            
            % Convert the box corners into the [x1 y1 x2 y2 x3 y3 x4 y4] 
            % format required by insertShape.
            bboxPolygon = reshape(bboxPoints', 1, []);            
            
            % Display a bounding box around the face being tracked.
            videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, 'LineWidth', 3);
            
            % Display tracked points.
            videoFrame = insertMarker(videoFrame, visiblePoints, '+', 'Color', 'white');
            
            % Reset the points.
            oldPoints = visiblePoints;
            setPoints(pointTracker, oldPoints);
        end
    end
    
    if ~isempty(bbox)
        x1 = bboxPoints(1,1);
        x2 = bboxPoints(2,1);
        x3 = bboxPoints(3,1);
        x4 = bboxPoints(4,1);
        y1 = bboxPoints(1,2);
        y2 = bboxPoints(2,2);
        y3 = bboxPoints(3,2);
        y4 = bboxPoints(4,2);
        dx1 = x2-x1;
        dx2 = x3-x4;
%         dy1 = y3-y2;
%         dy2 = y4-y1;
        x1 = x1 - 1 + 0.25*dx1;
        x2 = x2 - 0.25*dx1;
        x3 = x3 - 0.25*dx2;
        x4 = x4 + 0.25*dx2;
%         y3 = y3 - 0.8*dy1;
%         y4 = y4 - 0.8*dy2;
%         x=round([x1,x2,x3,x4]);
%         y=round([y1,y2,y3,y4]);
        x=[x1,x2,x3,x4];
        y=[y1,y2,y3,y4];
        x_double=double(x);
        y_double=double(y);
     mask=poly2mask(x_double,y_double,height,width); %finally working
     
    end
%     maskRed = redChannel.*mask;
%     maskGreen = greenChannel.*mask;
%     maskBlue = blueChannel.*mask;
    maskedImage = videoFrame;
    if ~isempty(bbox)
    for nChannel = 1:3
        colorChannel = maskedImage(:,:,nChannel);
        colorChannel = colorChannel.*mask;
        maskedImage(:,:,nChannel) = colorChannel;
    end
    end
    % Mask the image using bsxfun() function
    %maskedRgbImage = bsxfun(@times, rgbImage, cast(mask, 'like',
    %rgbImage)); %not yet try
    
    % Display the annotated video frame using the video player object.
    %videoPlayer(maskRed);
    %videoPlayer(maskGreen);
    %videoPlayer(maskBlue);
    videoPlayer(maskedImage);
    %videoPlayer(videoFrame);
%     
%     for nChannel = 1:3
%         colorChannel = maskedImage(:,:,nChannel);
%         colorChannel(colorChannel==0) = NaN;
%         maskedImage(:,:,nChannel) = colorChannel;
%     end
    for nChannel = 1:3
        colorChannel = maskedImage(:,:,nChannel);
        colorChannel(colorChannel==0) = NaN;
        rawColorSignal(nChannel,frameCount) =  mean(mean(colorChannel,'omitnan'),'omitnan');
    end
    % Mask ROI
%     maskRed(maskRed==0)=NaN;
%     maskGreen(maskGreen==0)=NaN;
%     maskBlue(maskBlue==0)=NaN;
    
    % Calculate mean
%     avgRed(frameCount) = mean(mean(maskRed,'omitnan'),'omitnan');
%     avgGreen(frameCount) = mean(mean(maskGreen,'omitnan'),'omitnan');
%     avgBlue(frameCount) = mean(mean(maskBlue,'omitnan'),'omitnan');
    
    % Check whether the video player window has been closed.
    runLoop = isOpen(videoPlayer);
    t=toc;
    t_tot=t_tot+t;
    t_hist(frameCount)=t;
    t_tot_hist(frameCount)=t_tot;
    trig = 11;
end

% Clean up.
clear cam;
release(videoPlayer);
release(pointTracker);
release(faceDetector);

%% References
% Viola, Paul A. and Jones, Michael J. "Rapid Object Detection using a
% Boosted Cascade of Simple Features", IEEE CVPR, 2001.
%
% Bruce D. Lucas and Takeo Kanade. An Iterative Image Registration 
% Technique with an Application to Stereo Vision. 
% International Joint Conference on Artificial Intelligence, 1981.
%
% Carlo Tomasi and Takeo Kanade. Detection and Tracking of Point Features. 
% Carnegie Mellon University Technical Report CMU-CS-91-132, 1991.
%
% Jianbo Shi and Carlo Tomasi. Good Features to Track. 
% IEEE Conference on Computer Vision and Pattern Recognition, 1994.
%
% Zdenek Kalal, Krystian Mikolajczyk and Jiri Matas. Forward-Backward
% Error: Automatic Detection of Tracking Failures.
% International Conference on Pattern Recognition, 2010
%samplingRate = round(totalFrame/t_tot);
samplingRate = 15;
lamda = 10;
%     if (samplingRate < 50) 
%       lambda = round(1.5*samplingRate)- 25;
%     else 
%       lambda = samplingRate;  % for 100 Hz
%     end
minFreq = 40;
maxFreq = 250;
dataLength = length(rawColorSignal);
heartRateBand = [minFreq maxFreq]*2/samplingRate;

  % detrend signals using Mean-Centering-And-Scaling technique [1]
  % (unconditionally since it's often a required step and it isn't harmful)
    w = samplingRate;    
    n = conv2(ones(3, dataLength), ones(1, w), 'same');
    meanIntensity = conv2(rawColorSignal, ones(1, w), 'same')./n;
    colorSignal = (rawColorSignal - meanIntensity)./meanIntensity;
    
   % additional detrending using smoothness prior approach [2]
   % we use crude approximation from the following table to compute lambda
    % SR:    | 30 | 50 | 100
    % lambda | 20 | 50 | 100
%     detrendingWindow = 10*samplingRate;
%     colorSignal = smoothness_priors_detrending(colorSignal, lambda, detrendingWindow);
  % use single green channel
    iPPG_1 = colorSignal(2, :);
  % use the difference of Green and Red   
    iPPG_2 = colorSignal(2, :) - colorSignal(1, :);
  % use CHROM method [1] 
    xs = 0.77*colorSignal(1, :) - 0.51*colorSignal(2, :);
    ys = 0.77*colorSignal(1, :) + 0.51*colorSignal(2, :) - 0.77*colorSignal(3, :);         
    alpha = std_ratio_sliding_win(xs, ys, fix(samplingRate));
    iPPG_3 = xs - alpha.*ys;
  % use POS method, see [a] for details and for the reference to the original paper
%     xs = colorSignal(2, :) - colorSignal(3, :);
%     ys = colorSignal(2, :) + colorSignal(3, :) - 2*colorSignal(1, :);  
%     alpha = std_ratio_sliding_win(xs, ys, fix(samplingRate));
%     iPPG_4 = xs + alpha.*ys; 
    
    % computes ratio of running standard deviations for two signals
function alpha = std_ratio_sliding_win(x, y, w)
  % n - number of elements in each window
  n = conv(ones(1, length(x)), ones(1, w), 'same');
  
  % compute running variation of x time series
  s = conv(x, ones(1, w), 'same');
  q = x.^ 2;    
  q = conv(q, ones(1, w), 'same');
  varXs = (q - s.^2 ./ n) ./ (n - 1);
    
  % compute running std of y time series  
  s = conv(y, ones(1, w), 'same');
  q = y.^ 2;    
  q = conv(q, ones(1, w), 'same');    
  varYs = (q - s.^2 ./ n)./(n - 1);  
    
  alpha = (varXs./varYs).^ 0.5; %compute std ratio
end