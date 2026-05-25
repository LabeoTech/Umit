function ImagesClassification(DataFolder, SaveFolder, BinningSpatial, BinningTemp, b_IgnoreStim, b_SubROI, chanName,trigPolarity, b_ApplyLPfilterToAnalogIn)
%IMAGESCLASSIFICATION Separate interleaved illumination frames into .dat files.
%
%   ImagesClassification(DataFolder, SaveFolder, BinningSpatial, ...
%       BinningTemp, b_IgnoreStim, b_SubROI, chanName, trigPolarity, ...
%       b_ApplyLPfilterToAnalogIn)
%
%   Separates raw Labeotech IOS img_*.bin files into one .dat/.mat pair per
%   output illumination channel. Raw frames are assumed to be stored in the
%   illumination sequence defined by info.txt. For example:
%
%       Red - Green - Amber - Fluo
%
%   is stored as:
%
%       R-G-A-F-R-G-A-F-R-G-A-F-...
%
%   Repeated illuminations are supported within each camera sequence. For
%   example:
%
%       Red - Green - Red
%
%   creates two output files:
%       red.dat   - contains both Red sequence positions in chronological order
%       green.dat - contains the Green sequence position
%
%   In this case, red.dat contains twice as many frames as green.dat before
%   temporal binning. The same grouping rule is applied independently to each
%   camera in dual-camera acquisitions. The channel frequency stored in the .mat metadata is
%   computed from the number of repetitions of each output channel:
%
%       Freq = FrameRateHz * NRepetitions / (NSequencePositions * BinningTemp)
%
%   The function uses the illumination sequence as the primary source of
%   truth. Each unique output channel is represented by one output group, and
%   each group stores the sequence positions assigned to that channel. Missing
%   frames are reconstructed in the raw stream, repeated sequence positions
%   are merged into one chronological output-channel stream, missing frames
%   are interpolated on that merged stream, and temporal binning is applied to
%   the final output-channel stream.
%
%   A ratio-aware finalization rule is used when writing files:
%
%       NFrames_channel / NRepetitions_channel
%
%   must be equal across all output channels from the same camera sequence.
%   Extra frames
%   belonging to incomplete final acquisition cycles are discarded before they
%   are written.
%
%   Repeated illumination channels are supported inside a single camera
%   sequence in MultiCam acquisitions. The unsupported case is a repeated
%   output channel across cameras, because the current .dat/.mat naming scheme
%   cannot represent the same channel tag independently for two cameras. In
%   that case, the function raises an error before writing channel files.
%
%   Inputs:
%       DataFolder                - Folder containing the raw Labeotech data.
%       SaveFolder                - Folder where output .dat/.mat files are saved.
%       BinningSpatial            - Spatial binning factor. Use 1 for no spatial
%                                   binning, 2 for 2x2 binning, 4 for 4x4, etc.
%       BinningTemp               - Temporal binning factor. Use 1 for no temporal
%                                   binning. Values > 1 merge consecutive frames
%                                   from each final output-channel stream.
%       b_IgnoreStim              - If true, stimulation signals are ignored.
%       b_SubROI                  - If true, the user selects or loads a spatial ROI.
%       chanName                  - Analog input channel(s) used for trigger/stim
%                                   extraction. Supported values include:
%                                       'Internal-main', 'Internal-Aux', 'AI1'...
%                                       'AI8', and 'StimDig'.
%       trigPolarity              - Trigger polarity. Use 'positive' or 'negative'.
%       b_ApplyLPfilterToAnalogIn - If true, low-pass filtering is applied to the
%                                   analog input during stimulation extraction.
%
%   Outputs:
%       One .dat file and one .mat metadata file per unique output channel.
%
%   Notes:
%       - Repeated channel positions are not averaged together unless temporal
%         binning explicitly combines them.
%       - For repeated channels, interpolation is performed after chronological
%         merging of all sequence positions assigned to the same output channel.
%       - The RAM usage remains chunk-based: img_*.bin files are memory-mapped
%         and processed one file at a time. Only small per-channel overflow
%         buffers are retained across chunks.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 6 || isempty(b_SubROI)
     b_SubROI = false;
end
if nargin < 7 || isempty(chanName)
     chanName = 'Internal-main';
end
if nargin < 8 || isempty(trigPolarity)
     trigPolarity = 'positive';
end
if nargin < 9 || isempty(b_ApplyLPfilterToAnalogIn)
     b_ApplyLPfilterToAnalogIn = false;
end
p = inputParser;
addRequired(p,'DataFolder',@(x)ischar(x)||isstring(x));
addRequired(p,'SaveFolder',@(x)ischar(x)||isstring(x));
addRequired(p,'BinningSpatial',@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addRequired(p,'BinningTemp',@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addRequired(p,'b_IgnoreStim',@(x)isscalar(x) && (isnumeric(x)||islogical(x)));
addOptional(p,'b_SubROI',false,@(x)isscalar(x) && (isnumeric(x)||islogical(x)));
addOptional(p,'chanName','',@(x) ( iscell(x) && ischar([x{:}]) ) || (ischar(x)||isstring(x)));
addOptional(p,'trigPolarity','positive',@(x)ischar(x)||isstring(x));
addOptional(p,'b_ApplyLPfilterToAnalogIn',false,@(x)isscalar(x) && (isnumeric(x)||islogical(x)));
parse(p,DataFolder,SaveFolder,BinningSpatial,BinningTemp,b_IgnoreStim,b_SubROI,chanName,trigPolarity,b_ApplyLPfilterToAnalogIn);

b_SubROI = p.Results.b_SubROI;
chanName = p.Results.chanName;
if ischar(chanName) || isstring(chanName)
    chanName = {chanName};
end
trigPolarity = char(p.Results.trigPolarity);
b_ApplyLPfilterToAnalogIn = p.Results.b_ApplyLPfilterToAnalogIn;


AcqInfoStream = ReadInfoFile(DataFolder);


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
            tmp = ReadAnalogsIn(DataFolder, SaveFolder, AcqInfoStream, chanName,trigPolarity,b_ApplyLPfilterToAnalogIn); %Read analog inputs.
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
% Check for each camera:
cam2_idx = contains({imgFilesList.name},'Cam2');
imgFileSet = {imgFilesList(~cam2_idx), imgFilesList(cam2_idx)};
imgFileSet(cellfun(@isempty,imgFileSet)) = [];
for ind = 1:length(imgFileSet)
    thisImgFilesList = imgFileSet{ind};
    % Check if all files exist:
    imgFileNames = sort({thisImgFilesList .name})';
    imgFileIndx = str2double(erase(imgFileNames, "img_" | "imgCam2_" | ".bin"));
    if imgFileIndx(1) ~=0  || any(diff(imgFileIndx) ~= 1)
        error('Image binary files (img_xxxxx.bin) missing! Classification aborted.')
    end
end

% Images files header description (see User Manual, page 26 for more
% details):
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
                load([pathname filesep filename]);%#ok
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
    
    LimX = [round(Pos(1)) round(Pos(1)+Pos(3)-1)];
    LimY = [round(Pos(2)) round(Pos(2)+Pos(4)-1)];
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
    Colors = struct('ID', {}, 'Color', {}, 'CamIdx', {}, 'FrameIdx', {}, 'Exposure', {}, 'Tag', {});
    for indC = 1:NbColors
        Colors(indC).ID = indC;
        eval(['Colors(' int2str(indC) ').Color = AcqInfoStream.Illumination' int2str(indC) '.Color;']);
        eval(['Colors(' int2str(indC) ').CamIdx = AcqInfoStream.Illumination' int2str(indC) '.CamIdx;']);
        eval(['Colors(' int2str(indC) ').FrameIdx = AcqInfoStream.Illumination' int2str(indC) '.FrameIdx;']);
        if( contains(Colors(indC).Color,'red', 'IgnoreCase', true) )
            Colors(indC).Tag = 'red';
            Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
        elseif( contains(Colors(indC).Color,'green', 'IgnoreCase', true) )
            Colors(indC).Tag = 'green';
            Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
        elseif( contains(Colors(indC).Color,'amber', 'IgnoreCase', true) )
            Colors(indC).Tag = 'yellow';
            Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
        elseif( contains(Colors(indC).Color,{'speckle'}, 'IgnoreCase', true) )
            Colors(indC).Tag = 'speckle';
            if( ~isfield(AcqInfoStream, 'ExposureSpeckleMsec') )
                Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
            else
                Colors(indC).Exposure = AcqInfoStream.ExposureSpeckleMsec;
            end
        else
            waveTag = regexp(Colors(indC).Color, '[0-9]{3}','match');
            if( ~isempty(waveTag) )
                Colors(indC).Tag = ['fluo_' waveTag{:}];
            else
                Colors(indC).Tag = 'fluo';
            end
            if( ~isfield(AcqInfoStream, 'ExposureFluoMsec') )
                Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
            else
                Colors(indC).Exposure = AcqInfoStream.ExposureFluoMsec;
            end
        end
    end

    % Repeated illuminations are supported within each camera sequence.
    % The repeated output channel across cameras is unsupported.
    
    CamList = unique([Colors.CamIdx]);
    TagList = {Colors.Tag};
    UniqueTags = unique(TagList, 'stable');
    for indT = 1:length(UniqueTags)
        idxTag = find(strcmp(TagList, UniqueTags{indT}));
        if( length(unique([Colors(idxTag).CamIdx])) > 1 )
            error(['Repeated illumination channel "' UniqueTags{indT} ...
                '" was detected across cameras. Repeated illuminations ' ...
                'are supported only within a single camera sequence.']);
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
    matHByCam = {matHcam1, matHcam2};
    matH = [matHcam1, matHcam2]; clear matHc*
    fprintf('Camera #2 Done. \n');
else
    Tags = fieldnames(AcqInfoStream);
    idx = contains(Tags, 'Illumination');
    NbColors = sum(idx);
    Colors = struct('ID', {}, 'Color', {}, 'Exposure', {}, 'Tag', {});
    for indC = 1:NbColors
        Colors(indC).ID = indC;
        eval(['Colors(' int2str(indC) ').Color = AcqInfoStream.Illumination' int2str(indC) '.Color;']);
        if( contains(Colors(indC).Color,'red', 'IgnoreCase', true) )
            Colors(indC).Tag = 'red';
            Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
        elseif( contains(Colors(indC).Color,'green', 'IgnoreCase', true) )
            Colors(indC).Tag = 'green';
            Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
        elseif( contains(Colors(indC).Color,'amber', 'IgnoreCase', true) )
            Colors(indC).Tag = 'yellow';
            Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
        elseif( contains(Colors(indC).Color,{'speckle'}, 'IgnoreCase', true) )
            Colors(indC).Tag = 'speckle';
            if( ~isfield(AcqInfoStream, 'ExposureSpeckleMsec') )
                Colors(indC).Exposure = AcqInfoStream.ExposureMsec;
            else
                Colors(indC).Exposure = AcqInfoStream.ExposureSpeckleMsec;
            end
        else
            waveTag = regexp(Colors(indC).Color, '[0-9]{3}','match');
            if( ~isempty(waveTag) )
                Colors(indC).Tag = ['fluo_' waveTag{:}];
            else
                Colors(indC).Tag = 'fluo';
            end
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

% Final output validation. ChannelsSort writes files using a ratio-aware
% rule, so this block only verifies the final output contract. In MultiCam
% mode, the contract is validated independently for each camera sequence.
if( AcqInfoStream.MultiCam )
    for indCam = 1:length(matHByCam)
        datLenList = cellfun(@(x) x.datLength,matHByCam{indCam});
        repeatCountList = cellfun(@(x) x.RepeatCount,matHByCam{indCam});
        baseLenList = datLenList ./ repeatCountList;
        assert(~(any(mod(baseLenList, 1) ~= 0) || ~isscalar(unique(baseLenList))), ...
            ['Channel classification failed final length validation for camera #' ...
            int2str(indCam) '. Expected datLength/RepeatCount to be equal ' ...
            'for all output channels from the same camera.']);
    end
else
    datLenList = cellfun(@(x) x.datLength,matH);
    repeatCountList = cellfun(@(x) x.RepeatCount,matH);
    baseLenList = datLenList ./ repeatCountList;
    assert(~(any(mod(baseLenList, 1) ~= 0) || ~isscalar(unique(baseLenList))), ...
        ['Channel classification failed final length validation. ' ...
        'Expected datLength/RepeatCount to be equal for all output channels.']);
end

% Remove internal grouping metadata from saved channel .mat files. These
% fields are only needed during classification and final validation.
cleanupVars = {'IlluminationSequenceIdx', 'RepeatCount'};
for ind = 1:length(matH)
    matFile = matH{ind}.Properties.Source;
    matVars = whos('-file', matFile);
    matVarNames = {matVars.name};
    varsToRemove = intersect(cleanupVars, matVarNames);

    if( ~isempty(varsToRemove) )
        matData = load(matFile);
        matData = rmfield(matData, varsToRemove);
        save(matFile, '-struct', 'matData');
    end
end

% --- LOCAL NESTED FUNCTION -----------------------------------------------

    function fColor = ChannelsSort(fList, colors)
        %CHANNELSSORT Split interleaved binary frames into grouped channel .dat files.
        %
        %   fColor = ChannelsSort(fList, colors) reads the image binaries in
        %   fList, separates frames using the illumination sequence in colors,
        %   merges repeated entries into chronological output-channel streams,
        %   applies interpolation/binning, and writes one .dat/.mat pair per
        %   unique output channel.
        %
        % Load Stim info:
        if( ~b_IgnoreStim )
            Stim = load([SaveFolder 'StimParameters.mat']);
        else
            Stim.Stim = 0;
        end
        stim_fn = fieldnames(Stim);
        stim_fn = stim_fn(startsWith(stim_fn, 'stim_', 'IgnoreCase', true));
        
        % Create one output group per unique output channel. Each group may
        % contain one or more illumination-sequence positions.
        fColor = {};
        fid = [];
        stimPos = 0;
        subNbColors = size(colors,2);
        OutputGroups = struct('Tag', {}, 'Color', {}, 'SeqIdx', {}, ...
            'RepeatCount', {}, 'Exposure', {}, 'hTag', {}, 'dTag', {});
        
        for indC = 1:subNbColors
            idxGroup = [];
            if( ~isempty(OutputGroups) )
                idxGroup = find(strcmp({OutputGroups.Tag}, colors(indC).Tag), 1, 'first');
            end
            if( isempty(idxGroup) )
                idxGroup = length(OutputGroups) + 1;
                OutputGroups(idxGroup).Tag = colors(indC).Tag;
                OutputGroups(idxGroup).Color = colors(indC).Color;
                OutputGroups(idxGroup).SeqIdx = indC;
                OutputGroups(idxGroup).RepeatCount = 1;
                OutputGroups(idxGroup).Exposure = colors(indC).Exposure;
                OutputGroups(idxGroup).hTag = [colors(indC).Tag '.mat'];
                OutputGroups(idxGroup).dTag = [colors(indC).Tag '.dat'];
            else
                if( colors(indC).Exposure ~= OutputGroups(idxGroup).Exposure )
                    error(['Repeated illumination channel "' colors(indC).Tag ...
                        '" has inconsistent exposure values.']);
                end
                OutputGroups(idxGroup).SeqIdx = [OutputGroups(idxGroup).SeqIdx indC];
                OutputGroups(idxGroup).RepeatCount = length(OutputGroups(idxGroup).SeqIdx);
            end
        end
        
        % Initialize output files and metadata.
        for indG = 1:length(OutputGroups)
            hTag = OutputGroups(indG).hTag;
            dTag = OutputGroups(indG).dTag;
            
            if( exist([SaveFolder hTag], 'file') )
                delete([SaveFolder hTag]);
            end
            fColor{indG} = matfile([SaveFolder hTag], 'Writable', true);
            fColor{indG}.datFile = dTag;
            fColor{indG}.datSize = [Ry, Rx]; % Flipped datSize
            fColor{indG}.Stim = Stim.Stim;
            for ii = 1:length(stim_fn)
                fColor{indG}.(stim_fn{ii}) = [];
            end
            fColor{indG}.datLength = 0;
            fColor{indG}.FirstDim = 'y';
            fColor{indG}.Datatype = 'single';
            fColor{indG}.datName = 'data';
            fColor{indG}.dim_names = {'Y', 'X', 'T'};
            fColor{indG}.Freq = AcqInfoStream.FrameRateHz * ...
                OutputGroups(indG).RepeatCount / (subNbColors * BinningTemp);
            fColor{indG}.tExposure = OutputGroups(indG).Exposure;
            fColor{indG}.IlluminationSequenceIdx = OutputGroups(indG).SeqIdx;
            fColor{indG}.RepeatCount = OutputGroups(indG).RepeatCount;
            fid(indG) = fopen([SaveFolder dTag],'w');
        end
        
        % Per-output buffers. ImBinBuffer stores unbinned temporal-overflow
        % frames. ImWriteBuffer stores binned frames that are waiting until
        % all channels can be committed while preserving datLength/RepeatCount.
        ImBinBuffer = cell(1,length(OutputGroups));
        StimBinBuffer = cell(1,length(OutputGroups));
        ImWriteBuffer = cell(1,length(OutputGroups));
        StimWriteBuffer = cell(1,length(OutputGroups));
        datLengthList = zeros(1,length(OutputGroups));
        for indG = 1:length(OutputGroups)
            ImBinBuffer{indG} = [];
            StimBinBuffer{indG} = [];
            ImWriteBuffer{indG} = [];
            StimWriteBuffer{indG} = [];
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
            
            % Keep only complete illumination cycles before assigning frames
            % to sequence positions. Partial cycles are carried to the next
            % binary file and discarded only if they remain after the last file.
            iData = cat(3, oIm, iData);
            overflow = mod(size(iData,3), subNbColors);
            if( overflow > 0 )
                oIm = iData(:,:, size(iData,3)-(overflow:-1:1)+1);
            else
                oIm = [];
            end
            Images = iData(:,:,1:(size(iData,3)-overflow));
            clear iData hData overflow;
            
            if isempty(Images)
                % No complete illumination cycle is available yet.
                fprintf('\n');
                continue
            end
            if( (~b_IgnoreStim) && (sum(Stim.Stim) ~= 0) && ~isempty(stim_fn) )
                SubStim = [];
                for ii = 1:length(stim_fn)
                    SubStim(ii,:) = Stim.(stim_fn{ii})(stimPos + (1:size(Images,3)));
                end
                stimPos = stimPos + size(SubStim,2);
            else
                SubStim = zeros(1,size(Images,3),'single');
            end
            
            Images = reshape(Images, ImRes_XY(2), ImRes_XY(1), subNbColors, []);
            SubStim = reshape(SubStim,size(SubStim,1), subNbColors, []);
            
            for indG = 1:length(OutputGroups)
                seqIdx = OutputGroups(indG).SeqIdx;
                
                % Merge all sequence positions assigned to this output group.
                % The reshape preserves chronological order within each cycle:
                % position 1, position 2, ..., then the next cycle.
                Ims = Images(:, :, seqIdx, :);
                Ims = reshape(Ims, ImRes_XY(2), ImRes_XY(1), []);
                SubStimTmp = SubStim(:, seqIdx, :);
                SubStimTmp = reshape(SubStimTmp, size(SubStim,1), []);
                
                %SubROI
                if( b_SubROI )
                    Ims = Ims(round(LimY(1)):round(LimY(2)),round(LimX(1)):round(LimX(2)),:);
                end
                
                if( ~isempty(ImBinBuffer{indG}) )
                    Ims = cat(3, ImBinBuffer{indG}, Ims);
                    SubStimTmp = [StimBinBuffer{indG}, SubStimTmp];
                end
                
                nFramesAll = size(Ims,3);
                nFullRaw = floor(nFramesAll/BinningTemp) * BinningTemp;
                
                if( nFullRaw < nFramesAll )
                    ImBinBuffer{indG} = Ims(:,:,nFullRaw+1:end);
                    StimBinBuffer{indG} = SubStimTmp(:,nFullRaw+1:end);
                else
                    ImBinBuffer{indG} = [];
                    StimBinBuffer{indG} = [];
                end
                
                if( nFullRaw == 0 )
                    continue
                end
                
                % Missing frames are interpolated after repeated positions have
                % been merged into the final output-channel chronology.
                imgSum = squeeze(sum(sum(Ims,1),2));
                imgSum = imgSum(:)';
                if( any(imgSum == 0) )
                    idx = find(imgSum > 1);
                    if( length(idx) > 1 )
                        nY = size(Ims,1);
                        nX = size(Ims,2);
                        Ims = interp1(idx, single(reshape(Ims(:,:,idx),[], length(idx)))', ...
                            1:size(Ims,3),'linear','extrap');
                        Ims = reshape(Ims', nY, nX, []);
                    elseif( isscalar(idx) )
                        Ims = repmat(single(Ims(:,:,idx)), 1, 1, size(Ims,3));
                    end
                end
                
                Ims = Ims(:,:,1:nFullRaw);
                SubStimTmp = SubStimTmp(:,1:nFullRaw);
                
                %Temporal Binning
                if( BinningTemp > 1 )
                    Ims = imresize3(Ims, [size(Ims,1), size(Ims,2),...
                        size(Ims,3)/BinningTemp], 'linear');
                    nStimRows = size(SubStimTmp,1);
                    SubStimTmp = reshape(SubStimTmp, nStimRows, BinningTemp, []);
                    SubStimTmp = reshape(ceil(mean(SubStimTmp, 2)), nStimRows, []);
                end
                
                %Spatial Binning
                if( BinningSpatial > 1 )
                    Ims = imresize(Ims,1/BinningSpatial);
                end
                
                if( isempty(ImWriteBuffer{indG}) )
                    ImWriteBuffer{indG} = single(Ims);
                    StimWriteBuffer{indG} = SubStimTmp;
                else
                    ImWriteBuffer{indG} = cat(3, ImWriteBuffer{indG}, single(Ims));
                    StimWriteBuffer{indG} = [StimWriteBuffer{indG}, SubStimTmp];
                end
            end
            clear Images SubStim;
            
            % Commit only frames that preserve the ratio-aware output length
            % invariant across all output groups.
            baseLenList = zeros(1,length(OutputGroups));
            for indG = 1:length(OutputGroups)
                if( isempty(ImWriteBuffer{indG}) )
                    nBufferedFrames = 0;
                else
                    nBufferedFrames = size(ImWriteBuffer{indG},3);
                end
                baseLenList(indG) = floor((datLengthList(indG) + nBufferedFrames) / ...
                    OutputGroups(indG).RepeatCount);
            end
            targetBaseLen = min(baseLenList);
            currentBaseLenList = datLengthList ./ [OutputGroups.RepeatCount];
            assert(~(any(mod(currentBaseLenList, 1) ~= 0) || ...
                    ~isscalar(unique(currentBaseLenList))), ...
                'Internal length bookkeeping failed during channel sorting.');
            currentBaseLen = currentBaseLenList(1);
            
            if( targetBaseLen > currentBaseLen )
                for indG = 1:length(OutputGroups)
                    nFramesToWrite = (targetBaseLen - currentBaseLen) * ...
                        OutputGroups(indG).RepeatCount;
                    if( nFramesToWrite > 0 )
                        fwrite(fid(indG), ImWriteBuffer{indG}(:,:,1:nFramesToWrite), 'single');
                        for ii = 1:length(stim_fn)
                            fColor{indG}.(stim_fn{ii}) = [fColor{indG}.(stim_fn{ii}); ...
                                StimWriteBuffer{indG}(ii,1:nFramesToWrite)'];
                        end
                        
                        if( nFramesToWrite < size(ImWriteBuffer{indG},3) )
                            ImWriteBuffer{indG} = ImWriteBuffer{indG}(:,:,nFramesToWrite+1:end);
                            StimWriteBuffer{indG} = StimWriteBuffer{indG}(:,nFramesToWrite+1:end);
                        else
                            ImWriteBuffer{indG} = [];
                            StimWriteBuffer{indG} = [];
                        end
                        datLengthList(indG) = datLengthList(indG) + nFramesToWrite;
                        fColor{indG}.datLength = datLengthList(indG);
                    end
                end
            end
            fprintf('\n');
        end
        
        % Discard incomplete final temporal-bin and ratio-overflow frames.
        for indG = 1:length(OutputGroups)
            if( ~isempty(ImBinBuffer{indG}) )
                fprintf('Discarding %d unbinned overflow frame(s) from %s.\n', ...
                    size(ImBinBuffer{indG},3), OutputGroups(indG).Tag);
            end
            if( ~isempty(ImWriteBuffer{indG}) )
                fprintf('Discarding %d ratio-overflow frame(s) from %s.\n', ...
                    size(ImWriteBuffer{indG},3), OutputGroups(indG).Tag);
            end
        end
        
        finalBaseLenList = datLengthList ./ [OutputGroups.RepeatCount];
        assert(~(any(mod(finalBaseLenList, 1) ~= 0) || ...
                ~isscalar(unique(finalBaseLenList))), ...
            ['Channel classification failed final length validation. ' ...
            'Expected datLength/RepeatCount to be equal for all output channels.']);
        
        for indG = 1:length(OutputGroups)
            fclose(fid(indG));
        end
    end
end
