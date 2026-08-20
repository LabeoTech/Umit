function exportBehavioralDataToVideo(RawFolder,filename)
%EXPORTBEHAVIORALDATATOVIDEO Reads camera frames from .bin behavioral camera and saves it to a video file.
%
%   This function reads the camera frames from .bin files in the specified 
%   raw folder and saves the frames to a video file. The video file format is 
%   determined by the extension of the filename.
%
%   Parameters:
%   RawFolder - A string specifying the path to the raw folder containing the .bin files.
%   filename - A string specifying the name of the output video file. 
%   The file format is determined by the extension of the filename. Supported formats are '.avi' and '.mp4'.
%
%   Example:
%   EXPORTBEHAVIORALDATATOVIDEO('C:\path\to\raw\folder', 'output.avi')
%
%   Note:
%   This function assumes that the behav_xxxxx.bin files and the info.txt
%   file are stored in the raw folder. 
%   
%   If the frame rate is higher than 160 Hz and the video format is '.mp4', 
%   the frame rate is reduced to 160 Hz to avoid errors.


binFileList = dir(fullfile(RawFolder,'behav_*.bin'));
if isempty(binFileList)
    warning('No behavioral camera data found in folder "%s"! Operation aborted!',RawFolder)
    return
end

% Get frame size from first file
fid = fopen(fullfile(binFileList(1).folder,binFileList(1).name),'r');
header = fread(fid, 5, 'int32'); fclose(fid);
imageSizeX = header(2);
imageSizeY = header(3);
frameFormat = {'int64',3,'frameHeader';'uint16',[imageSizeX imageSizeY],'cdata'};
% Initialize Video Writer object:
if endsWith(filename, '.avi')
    profile = 'Uncompressed AVI';
elseif endsWith(filename, '.mp4')
    profile = 'MPEG-4'; 
else
    error('Unsuported video file format');
end
writerObj = VideoWriter(filename, profile); % Name it.
% Set video parameters:
Info = ReadInfoFile(RawFolder);
writerObj.FrameRate = Info.BehavioralFrameRateHz;
if strcmpi(profile,'MPEG-4') && writerObj.FrameRate > 160
    warning('Frame rate reduced to 160Hz to avoid errors. To run at %dHz, save the file as .AVI format!',Info.BehavioralFrameRateHz)
    writerObj.FrameRate = 160;
end
open(writerObj);
for ii = 1:length(binFileList)
    fprintf('Reading file %s ...\n',binFileList(ii).name);
    fMap = memmapfile(fullfile(binFileList(ii).folder,binFileList(ii).name),...
        'Offset',4*5,'Format',frameFormat,'repeat',Inf);
    data = fMap.Data;
    fprintf('\tWriting data to video file...');
    writeVideo(writerObj,permute(double(cat(3,data.cdata))./double(intmax('uint16')),[2 1 4 3]));
    fprintf('Done.\n');
end
close(writerObj);
% Open the folder:
if ispc
    path = fileparts(filename);
    if isempty(path);path = pwd;end
    winopen(path);
end
disp('Done!');
                                                
end
