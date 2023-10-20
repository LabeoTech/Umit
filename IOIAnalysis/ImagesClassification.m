function ImagesClassification(DataFolder, SaveFolder, BinningSpatial, BinningTemp, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channels classification for Labeotech IOS systems.
%
% On IOS systems, acquisitions made with multiple channels (colors) are
% interlaced. For example, if the following channels were selected before
% recording: Red (R), Green (G), Yellow (Y) and Fluo (F)
% The saved images will be organised as follow (from frame #1):
% R-G-Y-F-R-G-Y-F-R-G-Y-F-R-G-Y-F-...-R-G-Y-F
% (see manual page 26 for more details)
%
% This function is used to separate each channel from an acquisition.
% At the end, 1 .dat file per channel will be generated.
% The file (.dat) contains the raw images in chronological order.
% The information about the acquisition (Freq., Stimulation vector, ROI, etc.)
% is stored in the "AcqInfos.mat" file in the SaveFolder.
%
%%% Input Parameters:
% 1- DataFolder Path:
%  Path contaning dataset from Labeo's system
% 2- SaveFolder Path:
%  Path where to save the .dat file(s).
% 3- Spatial Binning:
%  Set to 1 for no binning; 2 for a 2x2 binning; 4 for a 4x4 binning and so on...
% 4- Temporal Binning:
%  Set to 1 for no binning; Otherwise, enter de number of frames to be combined
%  for example: 4 will merge the images by group of 4. So, images #1,2,3,4
%  will become image #1 after avering them together. Images #5, 6, 7, 8 will
%  become frame #2, etc.
% Optional input parameters:
%  1- Region of Interest (ROI): This parameter is a boolean (0 or 1) to tell the software if
%       we want to keep the whole image or if we want to select a smaller ROI
%  2- backupFolder(char, default = ''): Name-value input pair.
%       Name of the subfolder where to move the data. Here, the user has 3
%       entry options:
%        - '' (empty string, default) : a dialog box will be presented to the User to ask what
%           to do with existing files.
%        - 'ERASE': All data from the SaveFolder are ERASED and no backup is done.
%        - 'AUTO': Lets the function create a subfolder with name "bkp_yyyymmddHHMMss"
%        - '<FOLDERNAME>': Type the name of the subfolder where to move the existing data before erasing
%            them from the SaveFolder.
%       See examples below:
%
%         % Example 1: Using the default backup folder value ('')
%         ImagesClassification(DataFolder,SaveFolder,BinningSpatial, BinningTemp,b_SubROI);
%         In this case, a dialog box will prompt the user to handle existing files.
%
%         % Example 2: Erasing existing data ('ERASE')
%         ImagesClassification(DataFolder,SaveFolder,BinningSpatial, BinningTemp,b_SubROI,'backupFolder','ERASE');
%         This will erase all existing data from the source folder without creating a backup.
%
%         % Example 3: Automatic backup folder creation
%         ImagesClassification(DataFolder,SaveFolder,BinningSpatial, BinningTemp,b_SubROI,'backupFolder','AUTO');
%         The function will create a subfolder with a name like "bkp_yyyymmddHHMMss" and move existing data there.
%
%         % Example 4: Specifying a custom backup folder
%         custom_folder = 'myBackupData';
%         ImagesClassification(DataFolder,SaveFolder,BinningSpatial, BinningTemp,b_SubROI,'backupFolder',custom_folder);
%         The function will move existing data to the specified custom backup folder.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Argument parsing:
p = inputParser;
addRequired(p,'DataFolder',@isfolder);
addRequired(p,'SaveFolder',@ischar);
addRequired(p,'BinningSpatial',@isscalar);
addRequired(p,'BinningTemp',@isscalar);
addOptional(p,'b_SubROI',false,@(x) islogical(x) || ismember(x,[0 1]));
addParameter(p,'backupFolder','',@ischar)
parse(p,DataFolder,SaveFolder,BinningSpatial,BinningTemp,varargin{:})
% Set Optional parameters:
b_SubROI = p.Results.b_SubROI;
backupFolder = p.Results.backupFolder;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Control for existing .dat files in the "SaveFolder":
% Excludes raw data, in case the SaveFolder is the same as the DataFolder.
fixedFiles = [dir(fullfile(SaveFolder,'*.bin'))... % binary files
    ;dir(fullfile(SaveFolder,'*nfo.txt'))... % info.txt and/or ExperimentInfo.txt
    ;dir(fullfile(SaveFolder,'Comments.txt'))... % Comments.txt
    ;dir(fullfile(SaveFolder,'Snapshot*.png'))]; % Snapshots
% Get list of files to move:
movFiles = dir(SaveFolder);
movFiles([movFiles.isdir] == 1) = [];
if ~isempty(fixedFiles)
    movFiles(ismember({movFiles.name},{fixedFiles.name})) = [];
end

if ~isempty(movFiles)
    
    % Create prompt if the moveToBackup optional parameter was not set:
    if isempty(backupFolder)
        
        % This function will erase all .dat and .mat files inside the
        % SaveFolder before the data import. If there are .dat files, the user will
        % be presented with the following options:
        %  - OPTION 1: Erase all data.
        %  - OPTION 2a: Move all folder content to a subfolder as backup (with the name "bkp_yyyymmddHHMMssFFF".
        %  - OPTION 2b: Move all folder content to a folder of the User's choice.
        
        choice = questdlg('The save folder already contains files. Please choose an option:', ...
            'Folder Contains Files', 'Erase all', 'Create backup', 'Cancel', 'Create backup');
        
        % Process the user's choice
        switch choice
            case 'Erase all'
                % Do not move to backup and erase all data from SaveFolder
                backupFolder = 'ERASE';
            case 'Create backup'
                % User chose to create a backup.
                answer = inputdlg('Type backup folder name:','BackupFolder',...
                    [1 60],{['bkp_' datestr(now(),'yyyymmddHHMMSS')]});
                if isempty(answer)
                    disp('Operation cancelled by User')
                    return
                elseif ~isempty(answer{:})
                    % Update backupfolder name
                    backupFolder = answer{:};
                else
                    % Force auto name to folder if the user provides an
                    % empty string:
                    backupFolder = ['bkp_' datestr(now(),'yyyymmddHHMMSS')];
                end
            otherwise
                % User chose to cancel the operation.
                disp('Operation cancelled by User')
                return
        end
        
    elseif strcmpi(backupFolder,'auto')
        % Auto name:
        backupFolder = ['bkp_' datestr(now(),'yyyymmddHHMMSS')];
    end
    
    if ~strcmpi(backupFolder,'erase')
        % Create subfolder
        if ~isfolder(fullfile(SaveFolder,backupFolder))
            mkdir(fullfile(SaveFolder,backupFolder));
        end
        % Copy data to subfolder
        arrayfun(@(x) copyfile(fullfile(x.folder,x.name), fullfile(x.folder,backupFolder,x.name)),movFiles)
    end
    % Delete all files from the current folder (exept .bin and .txt)
    arrayfun(@(x) delete(fullfile(x.folder,x.name)),movFiles);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AcqInfoStream = ReadInfoFile(DataFolder);

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = strcat(DataFolder, filesep);
end
if( ~strcmp(SaveFolder(end), filesep) )
    SaveFolder = strcat(SaveFolder, filesep);
end

if( ~isfield(AcqInfoStream, 'Camera_Model') )% For back compatibility with older versions
    AcqInfoStream.Camera_Model = 'D1024';     % Camera_Model was not used in former versions
end
% Create save folder if it does not exist.
if ~exist(SaveFolder,'dir')
    mkdir(SaveFolder)
end
% Change Spatial and temporal information in AcqInfoStream so the
% information refers to the output .dat files and not to the original .bin
% files.

% Update Binning information on AcqInfoStream:
AcqInfoStream.Binning = BinningSpatial;
% Update Sample rate in AcqInfoStream:
AcqInfoStream.FrameRateHz = AcqInfoStream.FrameRateHz/BinningTemp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Format and Header Information:

hWima = 5;
imgFilesList = dir([DataFolder 'img*.bin']);
%Images files header description (see User Manual, page 26 for more
%details):
header = memmapfile([DataFolder imgFilesList(1).name], ...
    'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
% Version = header.Data.header(1); %Data format version
nx = double(header.Data.header(2)); %Number of pixel along X axis
ny = double(header.Data.header(3)); %Number of pixel along Y axis
% FrameSz = header.Data.header(4); %Number of int32 saved for each image
% NbImsPefFile = single(header.Data.header(5)); %Number of images contained in each "img_" file.

%Header format for each individual image:
frameFormat = {'uint64', 3, 'framej';'uint16', [double(nx), double(ny)], 'imgj'};
ImRes_XY = [nx, ny];
SizeImage = nx*ny*2 + 3*8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SubROI...
if( b_SubROI )
    fprintf('Redefining Region Of Interest post-process: \n');
    %Dialog (there are different options to determine the new ROI):
    ButtonName = questdlg('Would you like to use a pre-defined ROI?', ...
        'ROI', ...
        'Pre-defined', 'Draw', 'Cancel', 'Draw');
    switch ButtonName  %Depending on user choice:
        case 'Pre-defined' %Used a ROI from an other acquisition:
            [filename, pathname] = uigetfile('*.mat', 'Select ROI file');
            if isequal(filename,0) || isequal(pathname,0)
                disp('User pressed cancel')
                Pos = [1 1 ImRes_XY(1) ImRes_XY(2)];
            else
                load([pathname filesep filename]);%#ok
            end
        case 'Draw' %Select ROI directly on a frame:
            dat = memmapfile([DataFolder...
                imgFilesList(1).name],...
                'Offset', hWima*4 + 5*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
            dat = rot90(fliplr(dat.Data.imgj));
            fig = figure('Name','Draw ROI','CloseRequestFcn',@closeFig); imagesc(dat); axis image;
            drawrectangle('Deletable',false,'Tag','myRectangle');
            title('Close figure to confirm')
            waitfor(fig)
            
        case 'Cancel' %User Changed is mind and want to use the original ROI
            disp('User pressed cancel')
            Pos = [1 1 ImRes_XY(1) ImRes_XY(2)];
    end
    
    LimX = [round(Pos(1)) round(Pos(1)+Pos(3)) - 1 ];
    LimY = [round(Pos(2)) round(Pos(2)+Pos(4)) - 1 ];
    save([SaveFolder 'ROI.mat'],'Pos'); %Save region of interest in a .mat file
    clear Pos
else
    LimX = [1 ImRes_XY(1)];
    LimY = [1 ImRes_XY(2)];
end
Rx = round((LimX(2) - LimX(1) + 1)/BinningSpatial);
Ry = round((LimY(2) - LimY(1) + 1)/BinningSpatial);
% Update Frame size in AcqInfoStream:
AcqInfoStream.Width = Rx;
AcqInfoStream.Height = Ry;
% Save AcqInfo:
save([SaveFolder 'AcqInfos.mat'],'AcqInfoStream');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How many colors and in which order?
fprintf('Sorting images per channels. \n');
if( AcqInfoStream.MultiCam )
    Tags = fieldnames(AcqInfoStream);
    idx = contains(Tags, 'Illumination');
    NbColors = sum(idx);
    Colors = struct('ID', {}, 'Color', {}, 'CamIdx', {}, 'FrameIdx', {}, 'Exposure', {});
    for indC = 1:NbColors
        Colors(indC).ID = indC;
        eval(['Colors(' int2str(indC) ').Color = AcqInfoStream.Illumination' int2str(indC) '.Color;']);
        eval(['Colors(' int2str(indC) ').CamIdx = AcqInfoStream.Illumination' int2str(indC) '.CamIdx;']);
        eval(['Colors(' int2str(indC) ').FrameIdx = AcqInfoStream.Illumination' int2str(indC) '.FrameIdx;']);
        if( contains(Colors(indC).Color,{'red', 'amber', 'green'}, 'IgnoreCase', true) )
            Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
        elseif( contains(Colors(indC).Color,{'speckle'}, 'IgnoreCase', true) )
            if( ~isfield(AcqInfoStream, 'ExposureSpeckleMsec') )
                Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
            else
                Colors(indC).Exposure = AcqInfoStream.ExposureSpeckleMsec;
            end
        else
            if( ~isfield(AcqInfoStream, 'ExposureFluoMsec') )
                Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
            else
                Colors(indC).Exposure = AcqInfoStream.ExposureFluoMsec;
            end
        end
    end
    
    %Camera 1:
    fprintf('Camera #1. \n');
    imgFilesList = dir([DataFolder 'img_*.bin']);
    idx = find(arrayfun(@(x) Colors(x).CamIdx == 1, 1:size(Colors,2)));
    [~, index] = sort([Colors(idx).FrameIdx]);
    idx = idx(index);
    datLenCam1 = ChannelsSort(imgFilesList, Colors(idx));
    fprintf('Camera #1 Done. \n');
    
    %Camera 2:
    fprintf('Camera #2. \n');
    imgFilesList = dir([DataFolder 'imgCam2_*.bin']);
    idx = find(arrayfun(@(x) Colors(x).CamIdx == 2, 1:size(Colors,2)));
    [~, index] = sort([Colors(idx).FrameIdx]);
    idx = idx(index);
    datLenCam2 = ChannelsSort(imgFilesList, Colors(idx));
    fprintf('Camera #2 Done. \n');
    datLen = vertcat(datLenCam1, datLenCam2); clear datLenCam*
else
    Tags = fieldnames(AcqInfoStream);
    idx = contains(Tags, 'Illumination');
    NbColors = sum(idx);
    Colors = struct('ID', {}, 'Color', {}, 'Exposure', {});
    for indC = 1:NbColors
        Colors(indC).ID = indC;
        eval(['Colors(' int2str(indC) ').Color = AcqInfoStream.Illumination' int2str(indC) '.Color;']);
        if( contains(Colors(indC).Color,{'red', 'amber', 'green'}, 'IgnoreCase', true) )
            Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
        elseif( contains(Colors(indC).Color,{'speckle'}, 'IgnoreCase', true) )
            if( ~isfield(AcqInfoStream, 'ExposureSpeckleMsec') )
                Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
            else
                Colors(indC).Exposure = AcqInfoStream.ExposureSpeckleMsec;
            end
        else
            if( ~isfield(AcqInfoStream, 'ExposureFluoMsec') )
                Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
            else
                Colors(indC).Exposure = AcqInfoStream.ExposureFluoMsec;
            end
        end
    end
    
    imgFilesList = dir([DataFolder 'img_*.bin']);
    datLen = ChannelsSort(imgFilesList, Colors);
    fprintf('Done. \n');
end

% Check if all color channels have the same length:
if ~isscalar(unique([datLen.Len]))
    disp('Fixing data length...')
    % Remove extra frames of channels so all have the same length
    newLen = min([datLen.Len]);
    datLen([datLen.Len] == newLen) = [];% Keep list of channels with extra frames.
    for ind = 1:length(datLen)
        fid = fopen(fullfile(SaveFolder,datLen(ind).File),'r');
        dat = fread(fid,Inf,'*single');fclose(fid);
        dat = reshape(dat,AcqInfoStream.Height,AcqInfoStream.Width,[]);
        dat = dat(:,:,1:newLen);
        fid = fopen(fullfile(SaveFolder,datLen(ind).File),'w');
        fwrite(fid,dat,'single');
        fclose(fid);
    end
    disp('Data length fixed.')
end


    function colorLen = ChannelsSort(fList, colors)
        % CHANNELSORT performs the classification of the raw data into the
        % existing colors and saves the data into separate .dat files.
        
        % For each color, initialise output files:
        
        fid = [];
        subNbColors = size(colors,2);
        colorLen = struct();
        for indC = 1:size(colors,2)
            
            if( contains(colors(indC).Color, {'red','green'},'IgnoreCase', true) )
                dTag = [lower(colors(indC).Color) '.dat'];
            elseif( contains(colors(indC).Color, 'amber', 'IgnoreCase', true) )
                dTag = 'yellow.dat';
            elseif( contains(colors(indC).Color, 'fluo', 'IgnoreCase', true) )
                waveTag = regexp(colors(indC).Color, '[0-9]{3}','match');
                if( ~isempty(waveTag) )
                    dTag = ['fluo_' waveTag{:} '.dat'];
                else
                    dTag = ['fluo' waveTag{:} '.dat'];
                end
            else
                dTag = 'speckle.dat';
            end
            fid(indC) = fopen([SaveFolder dTag],'w');
            %
            colorLen(indC).File = dTag;
            colorLen(indC).Len = 0;
        end
        
        % Opening Images Files:
        oIm = [];
        Cnt = 0;
        for indF = 1:size(fList,1)
            fprintf('Sorting %s.', fList(indF).name);
            data = memmapfile([DataFolder fList(indF).name],...
                'Offset', hWima*4, 'Format', frameFormat,...
                'repeat', inf);
            data = data.Data;
            hData = reshape([data.framej], 3, []);
            iData = reshape([data.imgj], ImRes_XY(1), ImRes_XY(2), []);
            iData = permute(iData,[2 1 3]);
            clear data;
            
            if( contains(AcqInfoStream.Camera_Model,{'D1024', 'D1312'}) )
                %                 if( indF == 1 )
                %                     hData = hData(:,(subNbColors + 1):end) - subNbColors;
                %                     iData = iData(:,:,(subNbColors + 1):end);
                %                 end
                SkipNFirst = sum(hData(1,:) == 0);
                MissingOffset = cumsum(hData(2,:));
                hData(1,:) = hData(1,:) + MissingOffset - hData(1,1) + 1;
                goodFrames = find(accumarray(hData(1, (SkipNFirst+1):end)',1)==1)';
                ConseqFromLeft = [1 diff(goodFrames,1,2)==1];
                ConseqFromRight = fliplr([true diff(fliplr(goodFrames),1,2)==-1]);
                goodFrames = goodFrames(ConseqFromLeft|ConseqFromRight);
                Images = zeros(ImRes_XY(2), ImRes_XY(1), (hData(1,end) - hData(1,1) + 1),'uint16');
                Images(:,:,goodFrames) = iData(:,:,goodFrames);
                iData = Images;
            elseif( contains(AcqInfoStream.Camera_Model, 'BFLY') )
                iNbF = hData(2,1) - Cnt;
                if( (hData(2,end) - hData(2,1)) > 0 )
                    fprintf('\t WARNING: %d missing frames.',(hData(2,end) - hData(2,1)));
                    Cnt = Cnt + (hData(2,end) - hData(2,1));
                end
                Images = zeros(ImRes_XY(2), ImRes_XY(1), (hData(1,end) - hData(1,1) + 1 + iNbF),'uint16');
                Images(:,:,(hData(1,:) - hData(1,1) + 1 + iNbF)) = iData;
                iData = Images;
                clear Images;
            elseif( any(hData(2,:)) ) %missing frames
                fprintf('\t WARNING: %d missing frames.',sum(hData(2,:)));
                hData(1,:) = 1:size(hData,2); %Ignore counter
                hData(1,:) = hData(1,:) + cumsum(hData(2,:));
                if( hData(2,1) >= 1 )
                    iNbF = hData(2,1);
                else
                    iNbF = 0;
                end
                Images = zeros(ImRes_XY(2), ImRes_XY(1), (hData(1,end) - hData(1,1) + 1 + iNbF),'uint16');
                
                Images(:,:,(hData(1,:) - hData(1,1) + 1 + iNbF)) = iData;
                iData = Images;
                clear Images;
            end
            iData = cat(3, oIm, iData);
            overflow = mod(size(iData,3), size(colors,2)*BinningTemp);
            if( overflow > 0 )
                oIm = iData(:,:, size(iData,3)-(overflow:-1:1)+1);
            else
                oIm = [];
            end
            Images = iData(:,:,1:(size(iData,3)-overflow));
            clear iData hData overflow;
            
            if isempty(Images)
                % In cases where the temporal binning causes the overflow
                % of all frames. This should happen only in the last ".bin"
                % file.
                break
            end
            Images = reshape(Images, ImRes_XY(2), ImRes_XY(1), subNbColors, []);
            
            for indC = 1:size(colors,2)
                Ims = squeeze(Images(:, :, indC, :));
                if( any(sum(sum(Ims,1),2) == 0) )
                    idx = find(sum(sum(Ims,1),2) > 1);
                    Ims = interp1(idx, single(reshape(Ims(:,:,idx),[], length(idx)))', 1:size(Ims,3),'linear','extrap');
                    Ims = reshape(Ims', ImRes_XY(2), ImRes_XY(1), []);
                end
                
                %SubROI
                if( b_SubROI )
                    Ims = Ims(round(LimY(1)):round(LimY(2)),round(LimX(1)):round(LimX(2)),:);
                end
                %Temporal Binning
                if( BinningTemp > 1 )
                    Ims = imresize3(Ims, [size(Ims,1), size(Ims,2),...
                        size(Ims,3)/BinningTemp], 'linear');
                end
                %Spatial Binning
                if( BinningSpatial > 1 )
                    Ims = imresize(Ims,1/BinningSpatial);
                end
                colorLen(indC).Len = colorLen(indC).Len + size(Ims,3);
                fwrite(fid(indC), single(Ims), 'single');
            end
            fprintf('\n');
        end
        for indC = 1:size(colors,2)
            fclose(fid(indC));
        end
    end

% Auxiliary functions:
    function closeFig(src,~)
        % CloseRequest function for Drawing ROI.
        % It saves the current position of the rectangle in the axis:
        rectH = findobj(src,'Tag','myRectangle');
        Pos = rectH.Position;
        delete(src)
    end

end

