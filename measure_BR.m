%% Track a Face in Scene
% Create System objects for reading and displaying video and for drawing a bounding 
% box of the object. 
%%clear all;

dataDir = './data';
vidFile = fullfile(dataDir,'PakIndra480.avi');
vid = VideoReader(vidFile);
len = vid.NumberOfFrames;
numberOfFrames = size(vid, 2);
fr = 0;
%% To obtain the time stamp where vision toolbox can't do
for frame = 1 : len
    video = read(vid,frame);
    fr= fr + 1;
    t(fr) = vid.CurrentTime;
end
clear vid;

videoFileReader = vision.VideoFileReader(vidFile);
%% Setting up configuration for color channel extraction
% numberOfFrames = size(videoFileReader, 2);
% meanRedLevels = zeros(numberOfFrames, 1);
% meanGreenLevels = zeros(numberOfFrames, 1);
% meanBlueLevels = zeros(numberOfFrames, 1);

%% 
% Read the first video frame, which contains the object, define the region.
%%
objectFrame = step(videoFileReader);
[height,width,~] = size(objectFrame);

%% 
% As an alternative, you can use the following commands to select the object 
% region using a mouse. The object must occupy the majority of the region:

figure(1); imshow(objectFrame);
pos=drawrectangle;
bbox=pos.Position;
% Show initial frame with a red bounding box.
%%
objectImage = insertShape(objectFrame,'Rectangle',bbox); 
figure(2);
imshow(objectImage);
title('Red box shows object region');
% Convert the first box into a list of 4 points
% This is needed to be able to visualize the rotation of the object.
bboxPoints = bbox2points(bbox(1, :));
%% 
% Detect interest points in the object region.
%%
points = detectMinEigenFeatures(rgb2gray(objectFrame),'ROI',bbox,'MinQuality',0.2);
%% 
% Display the detected points.
%%
pointImage = insertMarker(objectFrame,points.Location,'+','Color','white');
figure(3);
imshow(pointImage);
title('Detected interest points');
%% 
% Create a tracker object.
%%
tracker = vision.PointTracker('MaxBidirectionalError',2);
%% 
% Initialize the tracker with the initial point locations and the initial
% video frame.
points = points.Location;
%%
initialize(tracker, points, objectFrame);
%% 
% Read, track, display points, and results in each video frame.
videoPlayer  = vision.VideoPlayer('Position',...
    [100 100 [size(objectFrame, 2), size(objectFrame, 1)]+30]);
%%
% Make a copy of the points to be used for computing the geometric
% transformation between the points in the previous and the current frames
oldPoints = points;

%%
runLoop = true;
frameCount = 0;
[a,b]=size(points);
trace_pointY= zeros(a,len);
point_ref= points(:, 2);
rawColorSignal=zeros(3,frame);
reset(videoFileReader);
%%
while ~isDone(videoFileReader) && runLoop
      videoFrame = videoFileReader();
      [points, isFound] = tracker(videoFrame);
      frameCount = frameCount + 1;
      visiblePoints = points(isFound, :);
      oldInliers = oldPoints(isFound, :);
      trace_pointY(:, frameCount)= points(:, 2); %% for original position
      %trace_pointY(:, frameCount)= points(:, 1); %% For rotated (landscape) position (ex. resting hand)
      if size(visiblePoints, 1) >= 2 % need at least 2 points
        
        % Estimate the geometric transformation between the old points
        % and the new points and eliminate outliers
        [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
            oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);
        
        % Apply the transformation to the bounding box points
        bboxPoints = transformPointsForward(xform, bboxPoints);
 
        % Insert a bounding box around the object being tracked
        bboxPolygon = reshape(bboxPoints', 1, []);
        videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, ...
            'LineWidth', 3);
                
        % Display tracked points
        videoFrame = insertMarker(videoFrame, visiblePoints, '+', ...
            'Color', 'white');       
        
        % Reset the points
        oldPoints = visiblePoints;
        setPoints(tracker, oldPoints);
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
        dy1 = y3-y2;
        dy2 = y4-y1;
        x1 = x1 - 1 + 0.25*dx1;
        x2 = x2 - 0.25*dx1;
        x3 = x3 - 0.25*dx2;
        x4 = x4 + 0.25*dx2;
%         y3 = y3 - 0.8*dy1;
%         y4 = y4 - 0.8*dy2;
        x=round([x1,x2,x3,x4]);
        y=round([y1,y2,y3,y4]);
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
    for nChannel = 1:3
        colorChannel = maskedImage(:,:,nChannel);
        colorChannel(colorChannel==0) = NaN;
        rawColorSignal(nChannel,frameCount) =  mean(mean(colorChannel,'omitnan'),'omitnan');
    end
    % Display the annotated video frame using the video player object
    %frame=imresize(frame,0.5);
    videoPlayer(videoFrame);
    runLoop = isOpen(videoPlayer);

end

 %plot original
figure(4)
axes4 = axes(figure(4));
plot(t,(trace_pointY-point_ref));
xlim([0 60]);
xlabel('Time (second)');
ylabel('Pixel displacement');
set(axes4,'FontSize',24);
%plot(t,(trace_pointY-points(:, 1))) %plot rotated
%% 
% Release the video reader and player.
%%
release(videoPlayer);
release(videoFileReader);
release(tracker);
%% 
% Copyright 2012 The MathWorks, Inc.
samplingRate = 25;
    if (samplingRate < 50) 
      lambda = round(1.5*samplingRate)- 25;
    else 
      lambda = samplingRate;  % for 100 Hz
    end
minFreq = 40;
maxFreq = 250;
dataLength = length(rawColorSignal);
heartRateBand = [minFreq maxFreq]*2/samplingRate;

  % detrend signals using Mean-Centering-And-Scaling technique [1]
  % (unconditionally since it's often a required step and it isn't harmful)
    w = samplingRate;    
    n = conv2(ones(3, dataLength), ones(1, w), 'same');
    meanIntensity = conv2(rawColorSignal, ones(1, w), 'same')./n;
    colorSignal = double((rawColorSignal - meanIntensity)./meanIntensity);
    
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
    xs = colorSignal(2, :) - colorSignal(3, :);
    ys = colorSignal(2, :) + colorSignal(3, :) - 2*colorSignal(1, :);  
    alpha = std_ratio_sliding_win(xs, ys, fix(samplingRate));
    iPPG_4 = xs + alpha.*ys; 

%Normalize the PPG signal
% refinedPPG = movmean(refinedPPG, ippgSettings.maFilterLength);

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