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
%faceDetector = vision.CascadeObjectDetector('faceDetect.xml'); #own
%cascade classifier
faceDetector = vision.CascadeObjectDetector();
% Create the point tracker object.
pointTracker = vision.PointTracker('MaxBidirectionalError',0.5);

% Create the webcam object.
cam = imaq.VideoDevice('winvideo', 1, 'MJPG_640x480');
%cam = videoinput('winvideo',1, 'MJPG_640x480');
%cam.FramesPerTrigger = 1;
% Capture one frame to get its size.
videoFrame = step(cam);
frameSize = size(videoFrame);
[height,width,~] = size(videoFrame);
% Create the video player object. 
videoPlayer = vision.VideoPlayer('Position', [100 0 [frameSize(2), frameSize(1)]+30]);

%% Detection and Tracking
% Capture and process video frames from the webcam in a loop to detect and
% track a face. The loop will run for 400 frames or until the video player
% window is closed.

runLoop = true;
numPts = 0;
frameCount = 0;
t_tot = 0;

while runLoop && frameCount < 500
    tic;
    % Get the next frame.
    videoFrame = step(cam);
    videoFrameGray = rgb2gray(videoFrame);
    %[videoFrameGray,~,~] = rgb2hsv(videoFrame);
    frameCount = frameCount + 1;
    redChannel = (videoFrame(:, :, 1));
    greenChannel = (videoFrame(:, :, 2));
    blueChannel = (videoFrame(:, :, 3));
    
    if numPts < 15
        % Detection mode.
        bbox = faceDetector(videoFrameGray);
        
        if ~isempty(bbox)
            % Find corner points inside the detected region.
            points = detectMinEigenFeatures(videoFrameGray, 'ROI', bbox(1, :));
            
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
                oldInliers, visiblePoints, 'similarity', 'MaxDistance', 1.5);            
            
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
     x = bboxPoints(:,1);
     y = bboxPoints(:,2);

%      x_double=double(x');
%      y_double=double(y');
     
    x1 = bboxPoints(1,1);
    x2 = bboxPoints(2,1);
    x3 = bboxPoints(3,1);
    x4 = bboxPoints(4,1);
    y1 = bboxPoints(1,2);
    y2 = bboxPoints(2,2);
    y3 = bboxPoints(3,2);
    y4 = bboxPoints(4,2);
    
    fx = x1 - 1 + 0.3*(x2-x1);
    fy = y2 - 1 + 0.3*(y3-y2);
    lx = x2 - 0.3*(x2-x1);   %height includes first pixel
    ly = y3 - 0.3*(y3-y2);   %width includes first pixel
    
    x=[fx,lx,lx,fx];
    y=[fy,fy,ly,ly];
    
     x_double=double(x);
     y_double=double(y);
     
%     bbox_hist(:,frameCount) = bboxPolygon(1,:);
     mask=poly2mask(x_double,y_double,height,width); %finally working
     %mask(mask==0)=NaN;
    end
%     x1(frameCount) = bboxPolygon(:,1);
%     y1(frameCount) = bboxPolygon(:,2);
%     x2(frameCount) = bboxPolygon(:,3);
%     y2(frameCount) = bboxPolygon(:,4);
%     x3(frameCount) = bboxPolygon(:,5);
%     y3(frameCount) = bboxPolygon(:,6);
%     x4(frameCount) = bboxPolygon(:,7);
%     y4(frameCount) = bboxPolygon(:,8);
    
    %remember, bounding box is _x_ first and _x_ is column not row!
%     fr = round(bboxPolygon(:,2) + 0.3*(bboxPolygon(:,6)-bboxPolygon(:,2)));
%     fc = round(bboxPolygon(:,1) + 0.3*(bboxPolygon(:,3)-bboxPolygon(:,1)));
%     lr = round(bboxPolygon(:,6) - 0.3*(bboxPolygon(:,6)-bboxPolygon(:,2)));   %height includes first pixel
%     lc = round(bboxPolygon(:,3) - 0.3*(bboxPolygon(:,3)-bboxPolygon(:,1)));   %width includes first pixel
%     meanRedLevels(frameCount) = mean(mean(videoFrame(:, :, 1)));
%     meanGreenLevels(frameCount) = mean(mean(videoFrame(:, :, 2)));
%     meanBlueLevels(frameCount) = mean(mean(videoFrame(:, :, 3)));
%     meanRedLevels(frameCount) = mean(redChannel,[fr:lr; fc:lc]); %not working
%     meanGreenLevels(frameCount) = mean(greenChannel,[fr:lr; fc:lc]);
%     meanBlueLevels(frameCount) = mean(blueChannel,[fr:lr; fc:lc]); 

    maskRed = redChannel.*mask;
    maskGreen = greenChannel.*mask;
    maskBlue = blueChannel.*mask;
    
    % Display the annotated video frame using the video player object.
    %videoPlayer(maskRed);
    %videoPlayer(maskGreen);
    %videoPlayer(maskBlue);
    videoPlayer(videoFrame);
    
    maskRed(maskRed==0)=NaN;
    maskGreen(maskGreen==0)=NaN;
    maskBlue(maskBlue==0)=NaN;
    avgRed(frameCount) = mean(mean(maskRed,'omitnan'),'omitnan');
    avgGreen(frameCount) = mean(mean(maskGreen,'omitnan'),'omitnan');
    avgBlue(frameCount) = mean(mean(maskBlue,'omitnan'),'omitnan');
    % Check whether the video player window has been closed.
    runLoop = isOpen(videoPlayer);
    t=toc;
    t_tot=t_tot+t;
    t_hist(frameCount)=t;
    t_tot_hist(frameCount)=t_tot;
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
