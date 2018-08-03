function bins = parseVideo(video, thdblc, hdblc, tdblc, odblc, width, height)    
% Parse a video of the Franck lab counter.  Assumes that the camera is
% stable.
%
% Arguments are optional; if omitted, this function runs on the second video
% recorded.
%
% Parameters:
%
% video:  A string containing the name of the video file.  Must be a file
%         Matlab can decode.  I use avi files.
% thdblc: The coordinate of the bottom left corner of the thousands digit
%         in the video
% hdblc:  Same as above, but for the hundreds digit
% tdblc:  Same as above, but for the tens digit
% odblc:  Same as above, but for the ones digit
% width:  The width of the digit sub-image
% height: The height of the digit sub-image
%
% Happy parsing!
% 
% Example image reference creation line:
%
% MovieObj = VideoReader('long_run.avi');
% frame = read(MovieObj, 41);
% frame = rgb2gray(frame);
% ones(:,7) = reshape(frame(odblc(1):odblc(1)+height, odblc(2):odblc(2)+width), [18720 1]);
%
% Written by Andrew Parmet for Physics 3310, Spring 2016.

if (nargin == 0)
    video = 'long_run.avi';
    thdblc = [163 111];
    hdblc = [153 270];
    tdblc = [151 420];
    odblc = [143 594];
    width = 103;
    height = 179;
end


fprintf('Loading video; this may take a moment\n\n');

MovieObj = VideoReader(video);
numberOfFrames = get(MovieObj, 'NumberOfFrames');

% Matlab insists that you don't read frames after counting them, so load the
% video again...
MovieObj = VideoReader(video);

% col 1 is the time, col 2 is the count at that time
bins = zeros(3000, 2);

oneDimLength = (width+1)*(height+1);

ones = zeros(oneDimLength,10); % 1-9 and 0
tens = zeros(oneDimLength,10);
hundreds = zeros(oneDimLength,10); 
thousands = zeros(oneDimLength,10);

% Loads reference images; instantiates matrices ones, tens, hundreds, and
% thousands to contain column vectors where [digit](:,j) is the
% one-dimensional representation of the image of the digit j in the [digit]
% place on the counter. The image for 0 is in index 10.
load('basicData1.mat');

onesDistances = zeros(1,10);
tensDistances = zeros(1,10);
hundredsDistances = zeros(1,10);
thousandsDistances = zeros(1,10);

oneDimLength = (width+1)*(height+1);

currentBin = 1;

fprintf('Starting video parsing\n\n');

for k = 1:numberOfFrames
    
    frame = readFrame(MovieObj);
    frame = rgb2gray(frame);

    % ones digit check first; most likely to be bad.
    % Grab and reshape the subimages to be 1-dim vectors
    onesDigit = reshape(frame(odblc(1):odblc(1)+height, odblc(2):odblc(2)+width), [oneDimLength,1]);
    for i = 1:10
        onesDistances(i) = l2distance(double(onesDigit), double(ones(:,i)));
    end
        
    [minVal, idx] = min(onesDistances);
    sorted = sort(onesDistances);
    secondSmallest = sorted(2);
    
    % Image cleanliness condition
    if (1.5*minVal < secondSmallest)
        onesValue = mod(idx,10);
    else continue;
    end
        
    % tens digit check
    tensDigit = reshape(frame(tdblc(1):tdblc(1)+height, tdblc(2):tdblc(2)+width), [oneDimLength,1]);
    for i = 1:10
        tensDistances(i) = l2distance(double(tensDigit), double(tens(:,i)));
    end
        
    [minVal, idx] = min(tensDistances);
    sorted = sort(tensDistances);
    secondSmallest = sorted(2);
    
    % Image cleanliness condition
    if (1.5*minVal < secondSmallest)
        tensValue = mod(idx,10);
    else continue;
    end
    
    % hundreds digit check
    hundredsDigit = reshape(frame(hdblc(1):hdblc(1)+height, hdblc(2):hdblc(2)+width), [oneDimLength,1]);
    for i = 1:10
        hundredsDistances(i) = l2distance(double(hundredsDigit), double(hundreds(:,i)));
    end
        
    [minVal, idx] = min(hundredsDistances);
    sorted = sort(hundredsDistances);
    secondSmallest = sorted(2);
    
    % Image cleanliness condition
    if (1.5*minVal < secondSmallest)
        hundredsValue = mod(idx,10);
    else continue;
    end
    
    % thousands digit check
    thousandsDigit = reshape(frame(thdblc(1):thdblc(1)+height, thdblc(2):thdblc(2)+width), [oneDimLength,1]);  
    for i = 1:10
        thousandsDistances(i) = l2distance(double(thousandsDigit), double(thousands(:,i)));
    end
       
    [minVal, idx] = min(thousandsDistances);
    sorted = sort(thousandsDistances);
    secondSmallest = sorted(2);
    
    % Image cleanliness condition
    if (1.5*minVal < secondSmallest)
        thousandsValue = mod(idx,10);
    else continue;
    end
    
    % At this point, every digit is stable and we have a well-parsed
    % number.
    
    count = onesValue + 10*tensValue + 100*hundredsValue + 1000*thousandsValue;
    time = k/(30000/1001);  % compensate for frame rate
    
    bins(currentBin, 1) = time;
    bins(currentBin, 2) = count;
    
    currentBin = currentBin + 1;
    
    % Occassionally print status
    if (mod(currentBin, 50) == 0)
        fprintf('Current bin: %d\n', currentBin);
        fprintf('Current frame: %d of %d\n', k, numberOfFrames);
        fprintf('Current count: %d\n', count);
        fprintf('Current time: %3.1f\n\n', time);
    end
    
end
