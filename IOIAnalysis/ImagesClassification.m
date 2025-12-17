function ImagesClassification(DataFolder, SaveFolder, BinningSpatial, BinningTemp, b_IgnoreStim, b_SubROI, chanName,varargin)
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
% At the end, 1 .dat file and 1 .mat file per channel will be generated.
% The first file (.dat) contains the raw images in chronological order.
% The second file (.mat) contains all the informations about the
% acquisition (Freq., Stimulation vector, ROI, etc.)
%
%%% Input Parameters:
% 1- DataFolder Path:
%  Path contaning dataset from Labeo's system
% 2- SaveFolder Path:
%  Path where to save
% 3- Spatial Binning:
%  Set to 1 for no binning; 2 for a 2x2 binning; 4 for a 4x4 binning and so on...
% 4- Temporal Binning:
%  Set to 1 for no binning; Otherwise, enter de number of frames to be combined
%  for exemple: 4 will merge the images by group of 4. So, images #1,2,3,4
%  will become image #1 after avering them together. Images #5, 6, 7, 8 will
%  become frame #2, etc.
% 5- Ignore stimulation signal
%  boolean to tell the function if it should consider the
%  stimulation signal or not (0 = consider stim; 1= ignore stim)
% 6- Region of Interest (ROI)
%  this parameter is a boolean (0 or 1) to tell the software if
%  we want to keep the whole image or if we want to select a smaller ROI
% 7- Channel Name:
%   use this parameter to select a specific analog IN channel that contains
%   triggers. The available options are:
%       'Internal-main' (default): Main internal channel.
%       'Internal-Aux' : Auxiliary internal channel.
%       'AI1' : External Analog channel #1.
%       'AI2' : External Analog channel #2.
%       'AI3' : External Analog channel #3.
%       'AI4' : External Analog channel #4.
%       'AI5' : External Analog channel #5.
%       'AI6' : External Analog channel #6.
%       'AI7' : External Analog channel #7.
%       'AI8' : External Analog channel #8.
% 8- Trigger polarity (optional | default = positive): Set to "negative" if
%   the trigger onset and offset is marked by a falling and rising edges,
%   respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(varargin)
    trigPolarity = varargin{:};
else
    trigPolarity = 'positive';
end

if(nargin < 6)
    b_SubROI = 0;
end

AcqInfoStream = ReadInfoFile(DataFolder);


if strcmpi(chanName, 'Internal-Main')
    stimChan = 2;
elseif strcmpi(chanName, 'Internal-Aux')
    stimChan = 3;
else
    % Get position of the channel in the AnalogIN matrix:
    fn = fieldnames(AcqInfoStream);fn = fn(startsWith(fn,'AICh','IgnoreCase',true));
    if ~isempty(fn)
        % Look for the channel name directly from the info.txt file:
        chanNameList = cellfun(@(x) AcqInfoStream.(x),fn,'UniformOutput',false);
    else
        % Consider the following channel organization when there are no channel names in the info.txt file:
        chanNameList = {'CameraTrig','Internal-main', 'Internal-Aux','AI1', 'AI2','AI3','AI4','AI5','AI6','AI7','AI8'}; % List of existing Analog channel names.
    end
    [~,stimChan] = ismember(upper(chanName), upper(chanNameList));
    if stimChan == 0
        warning('Invalid channel name! The "Internal-main" channel will be read instead.');
        stimChan = 2;
    end
end

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = strcat(DataFolder, filesep);
end
if( ~strcmp(SaveFolder(end), filesep) )
    SaveFolder = strcat(SaveFolder, filesep);
end

if( ~isfield(AcqInfoStream, 'Camera_Model') ) %For back compatibility with older versions
    AcqInfoStream.Camera_Model = 'D1024';    %Camera_Model was not used in former versions
end
% Create save folder if it does not exist.
if ~exist(SaveFolder,'dir')
    mkdir(SaveFolder)
end
save([SaveFolder 'AcqInfos.mat'],'AcqInfoStream'); %To keep this information and avoid multiple reading of the txt file.

% Reading Stimulation Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulation parameters can be recovered from the analog inputs recordings
% when stimulation is generated from Labeo's acquisition system.

disp('Recovering stimulation parameters')
disp('**************************');

% tAIChan = AcqInfoStream.AINChannels; %Number of Analog Inputs

if( ~b_IgnoreStim ) %If user doesn't want to ignore Stimulation:
        
    if~isfield(AcqInfoStream, 'Stimulation') || ( AcqInfoStream.Stimulation == 0 ) 
        AcqInfoStream.Stimulation = 1;
    end
    Fields = fieldnames(AcqInfoStream); %Recovers Stimulation information from info.txt file
    idx = contains(Fields, 'stimulation','IgnoreCase',true);
    Fields = Fields(idx);
    
    for indS = 1:length(Fields)
        tmp = ReadAnalogsIn(DataFolder, SaveFolder, AcqInfoStream, stimChan,trigPolarity); %Read analog inputs.
        if(isnumeric(AcqInfoStream.(Fields{indS})) && ...
                AcqInfoStream.(Fields{indS}) > 0)
            break;
        end
    end
    if( isempty(tmp) ) %No stimulation found in info.txt file.
        fprintf('Stimulation not detected. \n');
        b_IgnoreStim = 1;
    end
    clear tmp
else
    fprintf('Stimulation ignored. \n');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Format and Header Information:

hWima = 5;
imgFilesList = dir([DataFolder 'img*.bin']); 
% Check if all files exist:
imgFileNames = sort({imgFilesList.name})';
imgFileIndx = str2double(erase(imgFileNames,"img_" | ".bin"));
if ~strcmpi(imgFileNames{1}, 'img_00000.bin') | any(diff(imgFileIndx)~=1)
    error('Image binary files missing! Classification aborted.')
end

% Images files header description (see User Manual, page 26 for more
% details):
header = memmapfile([DataFolder imgFilesList(1).name], ...
    'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);

Version = header.Data.header(1); %Data format version
nx = double(header.Data.header(2)); %Number of pixel along X axis
ny = double(header.Data.header(3)); %Number of pixel along Y axis
FrameSz = header.Data.header(4); %Number of int32 saved for each image
NbImsPefFile = single(header.Data.header(5)); %Number of images contained in each "img_" file.

%Header format for each individual image:
frameFormat = {'uint64', 3, 'framej';'uint16', [double(nx), double(ny)], 'imgj'};
ImRes_XY = [nx, ny];
SizeImage = nx*ny*2 + 3*8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SubROI...
if( b_SubROI )
    fprintf('Redifining Region Of Interest post-process: \n');
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
                load([pathname filesep filename]);
            end
        case 'Draw' %Select ROI directly on a frame:
            dat = memmapfile([DataFolder...
                imgFilesList(1).name],...
                'Offset', hWima*4 + 5*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
            dat = dat.Data.imgj;
            fig = figure; imagesc(dat);
            h = drawrectangle();
            wait(h);
            Pos = h.Position;
            close(fig);
        case 'Cancel' %User Changed is mind and want to use the original ROI
            disp('User pressed cancel')
            Pos = [1 1 ImRes_XY(1) ImRes_XY(2)];
    end
    
    LimX = [round(Pos(1)) round(Pos(1)+Pos(3))];
    LimY = [round(Pos(2)) round(Pos(2)+Pos(4))];
    save([SaveFolder 'ROI.mat'],'Pos'); %Save region of interest in a .mat file
else
    LimX = [1 ImRes_XY(1)];
    LimY = [1 ImRes_XY(2)];
end
Rx = round((LimX(2) - LimX(1) + 1)/BinningSpatial);
Ry = round((LimY(2) - LimY(1) + 1)/BinningSpatial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%How many colors and in which order?
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
    matHcam1 = ChannelsSort(imgFilesList, Colors(idx));
    fprintf('Camera #1 Done. \n');
    
    %Camera 2:
    fprintf('Camera #2. \n');
    imgFilesList = dir([DataFolder 'imgCam2_*.bin']);
    idx = find(arrayfun(@(x) Colors(x).CamIdx == 2, 1:size(Colors,2)));
    [~, index] = sort([Colors(idx).FrameIdx]);
    idx = idx(index);
    matHcam2 = ChannelsSort(imgFilesList, Colors(idx));
    matH = [matHcam1, matHcam2]; clear matHc*
    fprintf('Camera #2 Done. \n');
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
    matH = ChannelsSort(imgFilesList, Colors);
    fprintf('Done. \n');
end
% Check if all color channels have the same length:
datLenList = cellfun(@(x) x.datLength,matH);

if ~isscalar(unique(datLenList))
    disp('Fixing data length...')
    % Remove extra frames of channels so all have the same length
    newLen = min(datLenList);
    matH(datLenList == newLen) = [];% Keep list of channels with extra frames.
    for ind = 1:length(matH)
        fid = fopen(strrep(matH{ind}.Properties.Source,'.mat','.dat'),'r');
        dat = fread(fid,Inf,'*single');fclose(fid);
        dat = reshape(dat,matH{ind}.datSize(1,1),matH{ind}.datSize(1,2),[]);
        dat = dat(:,:,1:newLen);
        fid = fopen(strrep(matH{ind}.Properties.Source,'.mat','.dat'),'w');
        fwrite(fid,dat,'single');
        fclose(fid);
        matH{ind}.datLength = newLen;
    end
    
    disp('Data length fixed.')
end
    function fColor = ChannelsSort(fList, colors)
        % Load Stim info:
        if( ~b_IgnoreStim )
            Stim = load([SaveFolder 'StimParameters.mat']);
        else
            Stim.Stim = 0;
        end
        stim_fn = fieldnames(Stim);
        stim_fn = stim_fn(startsWith(stim_fn, 'stim_', 'IgnoreCase', true));
        
        % for each color, initialise output files:
        fColor = {};
        fid = [];
        stimPos = 0;
        subNbColors = size(colors,2);
        
        for indC = 1:size(colors,2)
            if( contains(colors(indC).Color, {'red','green'},'IgnoreCase', true) )
                hTag = [lower(colors(indC).Color) '.mat'];
                dTag = [lower(colors(indC).Color) '.dat'];
            elseif( contains(colors(indC).Color, 'amber', 'IgnoreCase', true) )
                hTag = ['yellow.mat'];
                dTag = ['yellow.dat'];
            elseif( contains(colors(indC).Color, 'fluo', 'IgnoreCase', true) )
                waveTag = regexp(colors(indC).Color, '[0-9]{3}','match');
                if( ~isempty(waveTag) )
                    hTag = ['fluo_' waveTag{:} '.mat'];
                    dTag = ['fluo_' waveTag{:} '.dat'];
                else
                    hTag = ['fluo' waveTag{:} '.mat'];
                    dTag = ['fluo' waveTag{:} '.dat'];
                end
            else
                hTag = 'speckle.mat';
                dTag = 'speckle.dat';
            end
            
            if( exist([SaveFolder hTag], 'file') )
                delete([SaveFolder hTag]);
            end
            fColor{indC} = matfile([SaveFolder hTag], 'Writable', true);
            fColor{indC}.datFile = dTag;
            fColor{indC}.datSize = [Ry, Rx]; % Flipped datSize
            fColor{indC}.Stim = Stim.Stim;
            for ii = 1:length(stim_fn)
                fColor{indC}.(stim_fn{ii}) = [];
            end
            fColor{indC}.datLength = 0;
            fColor{indC}.FirstDim = 'y';
            fColor{indC}.Datatype = 'single';
            fColor{indC}.datName = 'data';
            fColor{indC}.dim_names = {'Y', 'X', 'T'};
            fColor{indC}.Freq = (AcqInfoStream.FrameRateHz)/(size(colors,2)*BinningTemp);
            fColor{indC}.tExposure = colors(indC).Exposure;
            fid(indC) = fopen([SaveFolder dTag],'w');
        end
        
        %Opening Images Files:
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
                [~, goodFramesInIData] = ismember(goodFrames(1,:), hData(1,:));
                Images(:,:,goodFrames) = iData(:,:,goodFramesInIData);
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
            if( (~b_IgnoreStim) && (sum(Stim.Stim) ~= 0) )
                SubStim = [];
                for ii = 1:length(stim_fn)
                    SubStim(ii,:) = Stim.(stim_fn{ii})(stimPos + (1:size(Images,3)));
                end
                stimPos = stimPos + size(SubStim,2);
            else
                SubStim = zeros(1,size(Images,3),'single');
            end
            
            Images = reshape(Images, ImRes_XY(2), ImRes_XY(1), subNbColors, []);
            SubStim = reshape(SubStim,size(SubStim,1), subNbColors, BinningTemp, []);
            SubStim = ceil(mean(SubStim, 3));
            SubStim = reshape(SubStim, size(SubStim,1), subNbColors,[]);
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
                for ii = 1:length(stim_fn)
                    fColor{indC}.(stim_fn{ii}) = [fColor{indC}.(stim_fn{ii}); squeeze(SubStim(ii,indC,:))];
                end
                
                fwrite(fid(indC), single(Ims), 'single');
                fColor{indC}.datLength = fColor{indC}.datLength + size(Ims,3);
            end
            fprintf('\n');
        end
        for indC = 1:size(colors,2)
            fclose(fid(indC));
        end
    end
end
