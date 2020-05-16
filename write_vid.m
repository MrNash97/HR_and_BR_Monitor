    dataDir = './Data_PA2';
    vidFile = fullfile(dataDir,'ilham-dark.MOV');
    vid = VideoReader(vidFile);
    % Extract video info
    vidHeight = vid.Height;
    vidWidth = vid.Width;
    nChannels = 3;
    fr = vid.FrameRate;
    len = vid.NumberOfFrames;
    temp = struct('cdata', zeros(vidHeight, vidWidth, nChannels, 'uint8'), 'colormap', []);
    startIndex = 1;
    endIndex = len;
    
    vidOut = VideoWriter('ilham-dark.avi');
    vidOut.FrameRate = fr;

    open(vidOut)
    k = 0;
    for i=startIndex:endIndex
        k = k+1;
        temp.cdata = read(vid, i);
        [rgbframe,~] = frame2im(temp);
        rgbframe = imresize(rgbframe,[360 640]); %resize video
        frame = im2double(rgbframe);
        frame(frame > 1) = 1;
        frame(frame < 0) = 0;
        writeVideo(vidOut,im2uint8(frame));
    end

    disp('Finished')
    close(vidOut);

    