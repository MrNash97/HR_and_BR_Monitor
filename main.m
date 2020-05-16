function varargout = main(varargin)
% MAIN MATLAB code for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 16-Jun-2019 20:07:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)


% Choose default command line output for main
handles.output = hObject;
% Create video object
% Putting the object into manual trigger mode and then
% starting the object will make GETSNAPSHOT return faster
% since the connection to the camera will already have
% been established.
% 
% dataDir = './data';
% vidFile = fullfile(dataDir,'yudha-result.avi');
% handles.videoFileReader = vision.VideoFileReader(vidFile);
% handles.frame   = step(handles.videoFileReader);
% handles.fig=image(handles.frame);
rawColorSignal=[];
handles.rawColorSignal=rawColorSignal;
video = vision.VideoFileReader(); handles.video = video;
vidPlayer = vision.VideoPlayer; handles.vidPlayer = vidPlayer;
faceDetector = vision.CascadeObjectDetector(); handles.faceDetector = faceDetector;
pointTracker = vision.PointTracker('MaxBidirectionalError', 1); handles.pointTracker = pointTracker;
HsvMin = [0.00, 0.09, 0.34];handles.HsvMin = HsvMin;
HsvMax = [0.12, 0.52, 1.00];handles.HsvMax = HsvMax;
stdCoef = 1.5; handles.stdCoef = stdCoef;
numPts = 0; handles.numPts = numPts;
%frameCount = 0;
%handles.frameCount = frameCount;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles.output = hObject;
varargout{1} = handles.output;

% --- Executes on button press in Browse.
function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ video_file_name,video_file_path ] = uigetfile({'*.avi'},'Pick a video file');      %;*.png;*.yuv;*.bmp;*.tif'},'Pick a file');
if(video_file_path == 0)
    return;
end
input_video_file = [video_file_path,video_file_name];
fullpath = strcat(video_file_path,video_file_name);
set(handles.edit1,'String',fullpath);
axes(handles.axes1);set(handles.StartButton,'String','Start');
faceDetector=handles.faceDetector;
pointTracker=handles.pointTracker;
vid = VideoReader(input_video_file);
len = vid.NumberOfFrames;
clear vid;
video = vision.VideoFileReader(input_video_file);
%% load timestamp
h=info(video);
frameRate=h.VideoFrameRate;
rate=1/frameRate;
NumberOfFrames=len*rate;
time=rate;
vidPlayer = handles.vidPlayer;
vidFrame = step(video);
videoFrameGray = rgb2gray(vidFrame);
[height,width,~] = size(vidFrame);
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
            bboxPolygon = reshape(bboxPoints', 1, []);
            % Display a bounding box around the detected face.
            vidFrame = insertShape(vidFrame, 'Polygon', bboxPolygon, 'LineWidth', 3);
            % Display detected corners.
            vidFrame = insertMarker(vidFrame, xyPoints, '+', 'Color', 'white');
        end
if strcmp(get(handles.StartButton,'String'),'Start') && strcmp(get(handles.PlotButton,'String'),'Start')
    frameCount = 1;
    time = rate;
elseif strcmp(get(handles.StartButton,'String'),'Start')
    frameCount = 1;
    time = rate;
else
    frameCount = handles.frameCount;
    %time = handles.time;
end
step(vidPlayer,vidFrame);
imshow(videoFrameGray);
drawnow;
handles.frameCount = frameCount;
rawColorSignal=handles.rawColorSignal;
%%axis(handles.axes1,'off');
   for nChannel = 1:3
       colorChannel = vidFrame(:,:,nChannel);
       rawColorSignal(nChannel,frameCount) =  mean(mean(colorChannel));
   end
%plot(frameCount,rawColorSignal(1, :),frameCount,rawColorSignal(2, :),frameCount,rawColorSignal(3, :), handles.axes2);
%axes(handles.axes2)
plot(0:frameCount, rawColorSignal(1, :),'Parent',handles.axes2); drawnow;
%%axes(handles.axes1)
% Display Frame Number
%Update handles
handles.video = video;
handles.faceDetector = faceDetector;
handles.pointTracker = pointTracker;
handles.oldPoints = oldPoints;
handles.bboxPoints = bboxPoints;
handles.bbox = bbox;
handles.height = height;
handles.width = width;
handles.numPts = numPts;
handles.rate=rate;
handles.NumberOfFrames=NumberOfFrames;
handles.time=time;
guidata(hObject,handles);
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% --- Executes during object creation, after setting all properties.

function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in StartButton.
function StartButton_Callback(hObject, eventdata, handles)
% hObject    handle to StartButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.StartButton,'String'),'Pause')
    set(handles.StartButton,'String','Start');
else
    set(handles.StartButton,'String','Pause');
end
video = handles.video;
vidPlayer = handles.vidPlayer;
faceDetector=handles.faceDetector;
pointTracker=handles.pointTracker;
oldPoints=handles.oldPoints;
bboxPoints=handles.bboxPoints;
bbox=handles.bbox;
width=handles.width;
height=handles.height;
HsvMin=handles.HsvMin;
HsvMax=handles.HsvMax;
stdCoef=handles.stdCoef;
%NumberOfFrames=handles.NumberOfFrames;
rate=handles.rate;
if  isDone(video)
    reset(video)
    frameCount = 0;
    handles.frameCount = frameCount;
end
frameCount = handles.frameCount;
time = 0;
rawColorSignal=handles.rawColorSignal;
axes(handles.axes2)
 while ~isDone(video) && strcmp(get(handles.StartButton,'String'),'Pause')
    vidFrame   = step(video);
    step(vidPlayer,vidFrame);
    videoFrameGray = rgb2gray(vidFrame);
    frameCount = frameCount + 1;
    time = time + 0.02;
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
            vidFrame = insertShape(vidFrame, 'Polygon', bboxPolygon, 'LineWidth', 3);
            oldPoints = visiblePoints;
            setPoints(pointTracker, oldPoints);
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
    maskedImage = vidFrame;
    if ~isempty(bbox)
        maskedImage = maskedImage.*mask;
    end
    if ~isempty(bbox)
    maskedImage = maskedImage.*mask;
    % Add HSV filtering
    currentROI = ones(height,width);
    hsvFrame = rgb2hsv(maskedImage);
        for nChannel = 1:3
            hsvMask = roicolor(hsvFrame(:,:,nChannel), HsvMin(nChannel), HsvMax(nChannel));    
            currentROI = currentROI & hsvMask;
        end
        stdMask = currentROI;
        for nChannel = 1:3
            colorChannel = maskedImage(:,:,nChannel);
            channelMean = mean(mean(colorChannel(currentROI)));
            channelStd  = std2(colorChannel(currentROI));
            minChannelValue = channelMean - stdCoef*channelStd; 
            maxChannelValue = channelMean + stdCoef*channelStd;
            % we do not modify here currentROI since we need it to compute std corretly
            stdMask = stdMask & roicolor(colorChannel, minChannelValue, maxChannelValue);                
        end 
        currentROI = currentROI & stdMask;
        maskedImage = maskedImage.*currentROI;
    end
    imshow(maskedImage,'Parent',handles.axes3); %plot frame is specific axis
    drawnow;
    maskedImage(maskedImage==0) = NaN;
    for nChannel = 1:3
%        colorChannel = maskedImage(:,:,nChannel);
        rawColorSignal(nChannel,frameCount) =  mean(mean(maskedImage(:,:,nChannel),'omitnan'),'omitnan');
    end
    handles.frameCount = frameCount;
    handles.time = time;
    %Plotting data
    if strcmp(get(handles.PlotButton,'String'),'Stop') && handles.frameCount <= 45
        plot(1:handles.frameCount, rawColorSignal(1, :));
        %plot(0:handles.rate:time, rawColorSignal(1, :));
        xlim([2, 45]);
        drawnow;
        %grid on
    elseif strcmp(get(handles.PlotButton,'String'),'Stop') && handles.frameCount >= 45
        plot(1:handles.frameCount, rawColorSignal(1, :));
        xlim([(frameCount-43), (frameCount+5)]);
        drawnow;
        %grid on
    end
    
        %handles.time = time;
    %Plotting data
    %if strcmp(get(handles.PlotButton,'String'),'Stop') && handles.time <= 3
        %plot(1:handles.frameCount, rawColorSignal(1, :));
        %grid on
    %elseif strcmp(get(handles.PlotButton,'String'),'Stop') && handles.time >= 3
        %plot(0:handles.rate:handles.time, rawColorSignal(1, :));
        %xlim([(handles.time-2.96), (handles.time+0.1)]);
        %drawnow;
        %grid on
    %end
    
    vidFrame = insertMarker(vidFrame, visiblePoints, '+', 'Color', 'white');
    imshow(vidFrame,'Parent',handles.axes1); %plot frame is specific axis
    drawnow;
  end

%plot(frameCount,rawColorSignal(1, :),'r',frameCount,rawColorSignal(2, :),'g',frameCount,rawColorSignal(3, :),'b','Parent', handles.axes2);
%drawnow;
 set(handles.StartButton,'String','Start');

% --- Executes on button press in PlotButton.
function PlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.PlotButton,'String'),'Stop')
    set(handles.PlotButton,'String','Plot');
else
    set(handles.PlotButton,'String','Stop');
end
%surf(membrane(3))

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in fullplot.
function fullplot_Callback(hObject, eventdata, handles)
% hObject    handle to fullplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fullplot


% --- Executes on button press in moveplot.
function moveplot_Callback(hObject, eventdata, handles)
% hObject    handle to moveplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of moveplot
