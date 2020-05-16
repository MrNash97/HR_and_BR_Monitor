%% Face Detection and Tracking Using Live Video Acquisition
%% Setup
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
%videoFrame = step(cam);
videoFrame = snapshot(cam);
frameSize = size(videoFrame);
[height,width,~] = size(videoFrame);
% Create the video player object. 
videoPlayer = vision.VideoPlayer('Position', [100 0 [frameSize(2), frameSize(1)]+30]);
%% Detection and Tracking
trig = 0;
runLoop = true;
numPts = 0;
frameCount = 0;
t_tot = 0;
totalFrame = 500;
%t_hist=zeros(1,totalFrame);
t_tot_hist=zeros(1,totalFrame);
rawColorSignal=zeros(3,totalFrame);
% create circular array (overwriting value)
circularArraySz = 100;
t_hist          = zeros(1,circularArraySz); % create the circular cell array
nextIdx         = 1;                       % the index into the array in which
                                           % to insert the next matrix

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
            % Re-initialize the point tracker.
            xyPoints = points.Location;
            numPts = size(xyPoints,1);
            release(pointTracker);
            initialize(pointTracker, xyPoints, videoFrameGray);
            % Save a copy of the points.
            oldPoints = xyPoints;            
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
            bboxPolygon = reshape(bboxPoints', 1, []);            
            videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, 'LineWidth', 3);
            videoFrame = insertMarker(videoFrame, visiblePoints, '+', 'Color', 'white');
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
        x1 = x1 - 1 + 0.25*dx1;
        x2 = x2 - 0.25*dx1;
        x3 = x3 - 0.25*dx2;
        x4 = x4 + 0.25*dx2;
        x=[x1,x2,x3,x4];
        y=[y1,y2,y3,y4];
        x_double=double(x);
        y_double=double(y);
     mask=poly2mask(x_double,y_double,height,width); %finally working     
    end
    maskedImage = videoFrame;
    if ~isempty(bbox)
    for nChannel = 1:3
        colorChannel = maskedImage(:,:,nChannel);
        colorChannel = colorChannel.*mask;
        maskedImage(:,:,nChannel) = colorChannel;
    end
    end
    %videoPlayer(maskRed);
    %videoPlayer(maskGreen);
    %videoPlayer(maskBlue);
    videoPlayer(maskedImage);
    
    for nChannel = 1:3
        colorChannel = maskedImage(:,:,nChannel);
        colorChannel(colorChannel==0) = NaN;
        rawColorSignal(nChannel,frameCount) =  mean(mean(colorChannel,'omitnan'),'omitnan');
    end
    % Check whether the video player window has been closed.
    runLoop = isOpen(videoPlayer);
    t=toc;
    t_tot=t_tot+t;
    %t_hist(frameCount)=t;
    t_hist(nextIdx) = t;
    nextIdx = nextIdx + 1;
    t_tot_hist(frameCount)=t_tot;
    trig = 11;
    if nextIdx>circularArraySz
         nextIdx = 1;
     end
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