function varargout = ImagesClassification(DataFolder, SaveFolder, BinningSpatial, BinningTemp, varargin)
%IMAGESCLASSIFICATION Read interlaced binary acquisitions and split channels.
%
%   ImagesClassification(DataFolder, SaveFolder, BinningSpatial, BinningTemp)
%   ImagesClassification(DataFolder, SaveFolder, BinningSpatial, BinningTemp, b_SubROI)
%   ImagesClassification(..., 'backupOpts', backupOpts)
%   outFile = ImagesClassification(...)
%
%   This function reads raw interlaced .bin files produced by the LabeoTech
%   IOS systems, separates the illumination sequence into one .dat file per
%   unique output channel, and writes one shared AcqInfos.mat file in
%   SaveFolder.
%
%   Repeated illumination entries are supported within each camera sequence.
%   For example, Red-Green-Red creates red.dat and green.dat; red.dat stores
%   both Red sequence positions in chronological order and therefore has a
%   higher effective frame rate than green.dat. Repeated output channel tags
%   across cameras are rejected because the current file naming scheme cannot
%   represent the same output tag independently for multiple cameras.
%
%   Imported channel temporal metadata are stored in
%   AcqInfoStream.ImportedChannels. Spatial metadata remain at the top level
%   of AcqInfoStream.
%
%   Notes:
%       - The original function signature is intentionally preserved.
%       - Legacy per-channel metadata .mat files are no longer created.
%       - Existing files in SaveFolder are handled through genBackupFolder.
%       - If an output is requested, the function returns a file manifest in
%         varargout{1}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Argument parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'DataFolder', @isfolder);
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'BinningSpatial', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addRequired(p, 'BinningTemp', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addOptional(p, 'b_SubROI', false, @(x) islogical(x) || ismember(x, [0 1]));
addParameter(p, 'backupOpts', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));

parse(p, DataFolder, SaveFolder, BinningSpatial, BinningTemp, varargin{:});

DataFolder = char(string(p.Results.DataFolder));
SaveFolder = char(string(p.Results.SaveFolder));
BinningSpatial = p.Results.BinningSpatial;
BinningTemp = p.Results.BinningTemp;
b_SubROI = logical(p.Results.b_SubROI);
backupOpts = char(string(p.Results.backupOpts));

% Control for existing files and create backup or erase them.
genBackupFolder(SaveFolder, backupOpts,'eraseFolder',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outFile = {};
AcqInfoStream = ReadInfoFile(DataFolder);
rawFrameRateHz = AcqInfoStream.FrameRateHz;

if ~strcmp(DataFolder(end), filesep)
    DataFolder = strcat(DataFolder, filesep);
end
if ~strcmp(SaveFolder(end), filesep)
    SaveFolder = strcat(SaveFolder, filesep);
end

if ~isfield(AcqInfoStream, 'Camera_Model')
    AcqInfoStream.Camera_Model = 'D1024';
end
if ~exist(SaveFolder, 'dir')
    mkdir(SaveFolder)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Format and Header Information:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hWima = 5;
imgFilesList = dir([DataFolder 'img*.bin']);

% Check for each camera.
cam2_idx = contains({imgFilesList.name}, 'Cam2');
imgFileSet = {imgFilesList(~cam2_idx), imgFilesList(cam2_idx)};
imgFileSet(cellfun(@isempty, imgFileSet)) = [];
for ii = 1:length(imgFileSet)
    thisImgFilesList = imgFileSet{ii};
    imgFileNames = sort({thisImgFilesList.name})';
    imgFileIndx = str2double(erase(imgFileNames, "img_" | "imgCam2_" | ".bin"));
    if imgFileIndx(1) ~= 0 || any(diff(imgFileIndx) ~= 1)
        error('Image binary files (img_xxxxx.bin) missing! Classification aborted.')
    end
end

headerFiles = dir([DataFolder 'img_*.bin']);
if isempty(headerFiles)
    headerFiles = imgFilesList;
end
header = memmapfile([DataFolder headerFiles(1).name], ...
    'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
nx = double(header.Data.header(2));
ny = double(header.Data.header(3));
frameFormat = {'uint64', 3, 'framej'; 'uint16', [double(nx), double(ny)], 'imgj'};
ImRes_XY = [nx, ny];
SizeImage = nx * ny * 2 + 3 * 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SubROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if b_SubROI
    fprintf('Redefining Region Of Interest post-process:\n');
    ButtonName = questdlg('Would you like to use a pre-defined ROI?', ...
        'ROI', ...
        'Pre-defined', 'Draw', 'Cancel', 'Draw');
    switch ButtonName
        case 'Pre-defined'
            [filename, pathname] = uigetfile('*.mat', 'Select ROI file');
            if isequal(filename,0) || isequal(pathname,0)
                disp('User pressed cancel')
                Pos = [1 1 ImRes_XY(1) ImRes_XY(2)];
            else
                load([pathname filesep filename]); %#ok<LOAD>
            end
        case 'Draw'
            dat = memmapfile([DataFolder headerFiles(1).name], ...
                'Offset', hWima * 4 + 5 * SizeImage, ...
                'Format', frameFormat, 'repeat', 1);
            dat = rot90(fliplr(dat.Data.imgj));
            fig = figure('Name', 'Draw ROI', 'CloseRequestFcn', @closeFig);
            imagesc(dat);
            axis image;
            drawrectangle('Deletable', false, 'Tag', 'myRectangle');
            title('Close figure to confirm')
            waitfor(fig)
        case 'Cancel'
            disp('User pressed cancel')
            Pos = [1 1 ImRes_XY(1) ImRes_XY(2)];
    end

    LimX = [round(Pos(1)) round(Pos(1)+Pos(3)) - 1];
    LimY = [round(Pos(2)) round(Pos(2)+Pos(4)) - 1];
    save([SaveFolder 'ROI.mat'], 'Pos');
    clear Pos
else
    LimX = [1 ImRes_XY(1)];
    LimY = [1 ImRes_XY(2)];
end

Rx = round((LimX(2) - LimX(1) + 1) / BinningSpatial);
Ry = round((LimY(2) - LimY(1) + 1) / BinningSpatial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How many colors and in which order?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Sorting images per channels.\n');
if AcqInfoStream.MultiCam
    Tags = fieldnames(AcqInfoStream);
    idx = contains(Tags, 'Illumination');
    NbColors = sum(idx);
    Colors = struct('ID', {}, 'Color', {}, 'CamIdx', {}, 'FrameIdx', {}, 'Exposure', {}, 'Tag', {});
    for indC = 1:NbColors
        Colors(indC).ID = indC;
        eval(['Colors(' int2str(indC) ').Color = AcqInfoStream.Illumination' int2str(indC) '.Color;']);
        eval(['Colors(' int2str(indC) ').CamIdx = AcqInfoStream.Illumination' int2str(indC) '.CamIdx;']);
        eval(['Colors(' int2str(indC) ').FrameIdx = AcqInfoStream.Illumination' int2str(indC) '.FrameIdx;']);
        Colors(indC).Exposure = localResolveExposure(AcqInfoStream, Colors(indC).Color);
        Colors(indC).Tag = localGetChannelTag(Colors(indC).Color);
    end

    % Repeated output tags are allowed within one camera sequence but not
    % across cameras because each output tag maps to one .dat file name.
    TagList = {Colors.Tag};
    UniqueTags = unique(TagList, 'stable');
    for indT = 1:length(UniqueTags)
        idxTag = find(strcmp(TagList, UniqueTags{indT}));
        if length(unique([Colors(idxTag).CamIdx])) > 1
            error(['Repeated illumination channel "' UniqueTags{indT} ...
                '" was detected across cameras. Repeated illuminations ' ...
                'are supported only within a single camera sequence.']);
        end
    end

    fprintf('Camera #1.\n');
    imgFilesList = dir([DataFolder 'img_*.bin']);
    idx = find(arrayfun(@(x) Colors(x).CamIdx == 1, 1:size(Colors,2)));
    [~, index] = sort([Colors(idx).FrameIdx]);
    idx = idx(index);
    channelInfo1 = ChannelsSort(imgFilesList, Colors(idx));
    fprintf('Camera #1 Done.\n');

    fprintf('Camera #2.\n');
    imgFilesList = dir([DataFolder 'imgCam2_*.bin']);
    idx = find(arrayfun(@(x) Colors(x).CamIdx == 2, 1:size(Colors,2)));
    [~, index] = sort([Colors(idx).FrameIdx]);
    idx = idx(index);
    channelInfo2 = ChannelsSort(imgFilesList, Colors(idx));
    channelInfoByCam = {channelInfo1, channelInfo2};
    channelInfo = [channelInfo1, channelInfo2];
    fprintf('Camera #2 Done.\n');
else
    Tags = fieldnames(AcqInfoStream);
    idx = contains(Tags, 'Illumination');
    NbColors = sum(idx);
    Colors = struct('ID', {}, 'Color', {}, 'Exposure', {}, 'Tag', {});
    for indC = 1:NbColors
        Colors(indC).ID = indC;
        eval(['Colors(' int2str(indC) ').Color = AcqInfoStream.Illumination' int2str(indC) '.Color;']);
        Colors(indC).Exposure = localResolveExposure(AcqInfoStream, Colors(indC).Color);
        Colors(indC).Tag = localGetChannelTag(Colors(indC).Color);
    end

    imgFilesList = dir([DataFolder 'img_*.bin']);
    channelInfo = ChannelsSort(imgFilesList, Colors);
    channelInfoByCam = {channelInfo};
    fprintf('Done.\n');
end

% Final output validation. The output contract is validated independently
% for each camera sequence.
for indCam = 1:length(channelInfoByCam)
    datLenList = [channelInfoByCam{indCam}.datLength];
    repeatCountList = [channelInfoByCam{indCam}.RepeatCount];
    baseLenList = datLenList ./ repeatCountList;
    assert(~(any(mod(baseLenList, 1) ~= 0) || ~isscalar(unique(baseLenList))), ...
        ['Channel classification failed final length validation for camera #' ...
        int2str(indCam) '. Expected datLength/RepeatCount to be equal ' ...
        'for all output channels from the same camera.']);
end

% Save shared acquisition metadata for the processed outputs. Spatial
% properties stay at the top level. Imported-channel temporal metadata are
% stored in ImportedChannels.
baseLenList = [channelInfo.datLength] ./ [channelInfo.RepeatCount];
baseFreqList = [channelInfo.FrameRateHz] ./ [channelInfo.RepeatCount];
AcqInfoStream.Width = Rx;
AcqInfoStream.Height = Ry;
AcqInfoStream.Length = baseLenList(1);
AcqInfoStream.FrameRateHz = baseFreqList(1);
AcqInfoStream.BinningSpatial = BinningSpatial;
AcqInfoStream.BinningTemp = BinningTemp;
AcqInfoStream.Datatype = 'single';

if isfield(AcqInfoStream, 'ImportedChannels')
    AcqInfoStream = rmfield(AcqInfoStream, 'ImportedChannels');
end

for ind = 1:numel(channelInfo)
    importedInfo = struct();
    importedInfo.DatFile = channelInfo(ind).datFile;
    importedInfo.Tag = channelInfo(ind).Tag;
    importedInfo.Color = channelInfo(ind).Color;
    importedInfo.Length = channelInfo(ind).datLength;
    importedInfo.FrameRateHz = channelInfo(ind).FrameRateHz;
    importedInfo.ExposureMsec = channelInfo(ind).ExposureMsec;
    if ~isempty(channelInfo(ind).CamIdx)
        importedInfo.CamIdx = channelInfo(ind).CamIdx;
    end
    AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, importedInfo, 'Overwrite', true);
end

save([SaveFolder 'AcqInfos.mat'], 'AcqInfoStream');
outFile = [outFile, {channelInfo.datFile}];

if nargout ~= 0
    varargout{1} = outFile;
end

    function channelInfoOut = ChannelsSort(fList, colors)
        %CHANNELSSORT Split interleaved binary frames into grouped .dat files.

        subNbColors = size(colors, 2);
        OutputGroups = struct('Tag', {}, 'Color', {}, 'SeqIdx', {}, ...
            'RepeatCount', {}, 'Exposure', {}, 'dTag', {}, 'CamIdx', {});

        for indC = 1:subNbColors
            idxGroup = [];
            if ~isempty(OutputGroups)
                idxGroup = find(strcmp({OutputGroups.Tag}, colors(indC).Tag), 1, 'first');
            end

            if isempty(idxGroup)
                idxGroup = length(OutputGroups) + 1;
                OutputGroups(idxGroup).Tag = colors(indC).Tag;
                OutputGroups(idxGroup).Color = colors(indC).Color;
                OutputGroups(idxGroup).SeqIdx = indC;
                OutputGroups(idxGroup).RepeatCount = 1;
                OutputGroups(idxGroup).Exposure = colors(indC).Exposure;
                OutputGroups(idxGroup).dTag = [colors(indC).Tag '.dat'];
                if isfield(colors, 'CamIdx')
                    OutputGroups(idxGroup).CamIdx = colors(indC).CamIdx;
                else
                    OutputGroups(idxGroup).CamIdx = [];
                end
            else
                if colors(indC).Exposure ~= OutputGroups(idxGroup).Exposure
                    error(['Repeated illumination channel "' colors(indC).Tag ...
                        '" has inconsistent exposure values.']);
                end
                OutputGroups(idxGroup).SeqIdx = [OutputGroups(idxGroup).SeqIdx indC];
                OutputGroups(idxGroup).RepeatCount = length(OutputGroups(idxGroup).SeqIdx);
            end
        end

        channelInfoOut = repmat(struct( ...
            'datFile', '', ...
            'datSize', [Ry, Rx], ...
            'datLength', 0, ...
            'Tag', '', ...
            'Color', '', ...
            'FrameRateHz', [], ...
            'ExposureMsec', [], ...
            'CamIdx', [], ...
            'RepeatCount', []), ...
            1, length(OutputGroups));
        fid = zeros(1, length(OutputGroups));

        for indG = 1:length(OutputGroups)
            dTag = OutputGroups(indG).dTag;
            channelInfoOut(indG).datFile = dTag;
            channelInfoOut(indG).datSize = [Ry, Rx];
            channelInfoOut(indG).datLength = 0;
            channelInfoOut(indG).Tag = OutputGroups(indG).Tag;
            channelInfoOut(indG).Color = OutputGroups(indG).Color;
            channelInfoOut(indG).FrameRateHz = rawFrameRateHz * ...
                OutputGroups(indG).RepeatCount / (subNbColors * BinningTemp);
            channelInfoOut(indG).ExposureMsec = OutputGroups(indG).Exposure;
            channelInfoOut(indG).CamIdx = OutputGroups(indG).CamIdx;
            channelInfoOut(indG).RepeatCount = OutputGroups(indG).RepeatCount;
            if exist([SaveFolder dTag], 'file')
                delete([SaveFolder dTag]);
            end
            fid(indG) = fopen([SaveFolder dTag], 'w');
        end

        ImBinBuffer = cell(1, length(OutputGroups));
        ImWriteBuffer = cell(1, length(OutputGroups));
        datLengthList = zeros(1, length(OutputGroups));
        for indG = 1:length(OutputGroups)
            ImBinBuffer{indG} = [];
            ImWriteBuffer{indG} = [];
        end

        oIm = [];
        Cnt = 0;
        for indF = 1:size(fList,1)
            fprintf('Sorting %s.', fList(indF).name);
            data = memmapfile([DataFolder fList(indF).name], ...
                'Offset', hWima * 4, 'Format', frameFormat, 'repeat', inf);
            data = data.Data;
            hData = reshape([data.framej], 3, []);
            iData = reshape([data.imgj], ImRes_XY(1), ImRes_XY(2), []);
            iData = permute(iData, [2 1 3]);
            clear data;

            if contains(AcqInfoStream.Camera_Model, {'D1024', 'D1312'})
                SkipNFirst = sum(hData(1,:) == 0);
                MissingOffset = cumsum(hData(2,:));
                hData(1,:) = hData(1,:) + MissingOffset - hData(1,1) + 1;
                goodFrames = find(accumarray(hData(1, (SkipNFirst+1):end)',1) == 1)';
                ConseqFromLeft = [1 diff(goodFrames,1,2) == 1];
                ConseqFromRight = fliplr([true diff(fliplr(goodFrames),1,2) == -1]);
                goodFrames = goodFrames(ConseqFromLeft | ConseqFromRight);
                Images = zeros(ImRes_XY(2), ImRes_XY(1), (hData(1,end) - hData(1,1) + 1), 'uint16');
                [~, goodFramesInIData] = ismember(goodFrames(1,:), hData(1,:));
                Images(:,:,goodFrames) = iData(:,:,goodFramesInIData);
                iData = Images;
            elseif contains(AcqInfoStream.Camera_Model, 'BFLY')
                iNbF = hData(2,1) - Cnt;
                if (hData(2,end) - hData(2,1)) > 0
                    fprintf('\t WARNING: %d missing frames.', (hData(2,end) - hData(2,1)));
                    Cnt = Cnt + (hData(2,end) - hData(2,1));
                end
                Images = zeros(ImRes_XY(2), ImRes_XY(1), (hData(1,end) - hData(1,1) + 1 + iNbF), 'uint16');
                Images(:,:,(hData(1,:) - hData(1,1) + 1 + iNbF)) = iData;
                iData = Images;
                clear Images;
            elseif any(hData(2,:))
                fprintf('\t WARNING: %d missing frames.', sum(hData(2,:)));
                hData(1,:) = 1:size(hData,2);
                hData(1,:) = hData(1,:) + cumsum(hData(2,:));
                if hData(2,1) >= 1
                    iNbF = hData(2,1);
                else
                    iNbF = 0;
                end
                Images = zeros(ImRes_XY(2), ImRes_XY(1), (hData(1,end) - hData(1,1) + 1 + iNbF), 'uint16');
                Images(:,:,(hData(1,:) - hData(1,1) + 1 + iNbF)) = iData;
                iData = Images;
                clear Images;
            end

            iData = cat(3, oIm, iData);
            overflow = mod(size(iData,3), subNbColors);
            if overflow > 0
                oIm = iData(:,:, size(iData,3) - (overflow:-1:1) + 1);
            else
                oIm = [];
            end
            Images = iData(:,:,1:(size(iData,3) - overflow));
            clear iData hData overflow;

            if isempty(Images)
                fprintf('\n');
                continue
            end
            Images = reshape(Images, ImRes_XY(2), ImRes_XY(1), subNbColors, []);

            for indG = 1:length(OutputGroups)
                seqIdx = OutputGroups(indG).SeqIdx;
                Ims = Images(:, :, seqIdx, :);
                Ims = reshape(Ims, ImRes_XY(2), ImRes_XY(1), []);

                if b_SubROI
                    Ims = Ims(round(LimY(1)):round(LimY(2)), round(LimX(1)):round(LimX(2)), :);
                end

                if ~isempty(ImBinBuffer{indG})
                    Ims = cat(3, ImBinBuffer{indG}, Ims);
                end

                nFramesAll = size(Ims,3);
                nFullRaw = floor(nFramesAll / BinningTemp) * BinningTemp;

                if nFullRaw < nFramesAll
                    ImBinBuffer{indG} = Ims(:,:,nFullRaw+1:end);
                else
                    ImBinBuffer{indG} = [];
                end

                if nFullRaw == 0
                    continue
                end

                imgSum = squeeze(sum(sum(Ims,1),2));
                imgSum = imgSum(:)';
                if any(imgSum == 0)
                    idxGood = find(imgSum > 1);
                    if length(idxGood) > 1
                        nY = size(Ims,1);
                        nX = size(Ims,2);
                        Ims = interp1(idxGood, single(reshape(Ims(:,:,idxGood),[], length(idxGood)))', ...
                            1:size(Ims,3), 'linear', 'extrap');
                        Ims = reshape(Ims', nY, nX, []);
                    elseif length(idxGood) == 1
                        Ims = repmat(single(Ims(:,:,idxGood)), 1, 1, size(Ims,3));
                    end
                end

                Ims = Ims(:,:,1:nFullRaw);

                if BinningTemp > 1
                    Ims = imresize3(Ims, [size(Ims,1), size(Ims,2), ...
                        size(Ims,3) / BinningTemp], 'linear');
                end

                if BinningSpatial > 1
                    Ims = imresize(Ims, 1 / BinningSpatial);
                end

                if isempty(ImWriteBuffer{indG})
                    ImWriteBuffer{indG} = single(Ims);
                else
                    ImWriteBuffer{indG} = cat(3, ImWriteBuffer{indG}, single(Ims));
                end
            end
            clear Images;

            baseLenList = zeros(1, length(OutputGroups));
            for indG = 1:length(OutputGroups)
                if isempty(ImWriteBuffer{indG})
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

            if targetBaseLen > currentBaseLen
                for indG = 1:length(OutputGroups)
                    nFramesToWrite = (targetBaseLen - currentBaseLen) * ...
                        OutputGroups(indG).RepeatCount;
                    if nFramesToWrite > 0
                        fwrite(fid(indG), ImWriteBuffer{indG}(:,:,1:nFramesToWrite), 'single');

                        if nFramesToWrite < size(ImWriteBuffer{indG},3)
                            ImWriteBuffer{indG} = ImWriteBuffer{indG}(:,:,nFramesToWrite+1:end);
                        else
                            ImWriteBuffer{indG} = [];
                        end
                        datLengthList(indG) = datLengthList(indG) + nFramesToWrite;
                        channelInfoOut(indG).datLength = datLengthList(indG);
                    end
                end
            end
            fprintf('\n');
        end

        for indG = 1:length(OutputGroups)
            if ~isempty(ImBinBuffer{indG})
                fprintf('Discarding %d unbinned overflow frame(s) from %s.\n', ...
                    size(ImBinBuffer{indG},3), OutputGroups(indG).Tag);
            end
            if ~isempty(ImWriteBuffer{indG})
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

    function closeFig(src,~)
        % CloseRequest function for drawing ROI.
        rectH = findobj(src, 'Tag', 'myRectangle');
        Pos = rectH.Position;
        delete(src)
    end
end

function exposure = localResolveExposure(AcqInfoStream, colorName)
%LOCALRESOLVEEXPOSURE Resolve per-channel exposure metadata.

if contains(colorName, {'red','amber','green'}, 'IgnoreCase', true)
    exposure = AcqInfoStream.ExposureMsec;
elseif contains(colorName, {'speckle'}, 'IgnoreCase', true)
    if ~isfield(AcqInfoStream, 'ExposureSpeckleMsec')
        exposure = AcqInfoStream.ExposureMsec;
    else
        exposure = AcqInfoStream.ExposureSpeckleMsec;
    end
else
    if ~isfield(AcqInfoStream, 'ExposureFluoMsec')
        exposure = AcqInfoStream.ExposureMsec;
    else
        exposure = AcqInfoStream.ExposureFluoMsec;
    end
end
end

function tag = localGetChannelTag(colorName)
%LOCALGETCHANNELTAG Convert illumination color metadata to output tag.

if contains(colorName, 'red', 'IgnoreCase', true)
    tag = 'red';
elseif contains(colorName, 'green', 'IgnoreCase', true)
    tag = 'green';
elseif contains(colorName, 'amber', 'IgnoreCase', true)
    tag = 'yellow';
elseif contains(colorName, 'speckle', 'IgnoreCase', true)
    tag = 'speckle';
else
    waveTag = regexp(colorName, '[0-9]{3}', 'match');
    if ~isempty(waveTag)
        tag = ['fluo_' waveTag{:}];
    else
        tag = 'fluo';
    end
end
end
