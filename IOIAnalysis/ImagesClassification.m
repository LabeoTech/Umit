function varargout = ImagesClassification(DataFolder, SaveFolder, BinningSpatial, BinningTemp, b_SubROI)
%IMAGESCLASSIFICATION Read interlaced binary acquisitions and split channels.
%
%   ImagesClassification(DataFolder, SaveFolder, BinningSpatial, BinningTemp, b_SubROI)
%   outFile = ImagesClassification(...)
%
%   This function reads raw interlaced .bin files produced by the LabeoTech
%   IOS systems, separates the channels, writes one .dat file per channel,
%   and writes one shared AcqInfos.mat file in SaveFolder.
%
%   Notes:
%       - The original function signature is intentionally preserved.
%       - Legacy per-channel metadata .mat files are no longer created.
%       - If an output is requested, the function returns a file manifest in
%         varargout{1}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    b_SubROI = 0;
end

outFile = {};
AcqInfoStream = ReadInfoFile(DataFolder);

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

hWima = 5;
imgFilesList = dir([DataFolder 'img*.bin']);
imgFileNames = sort({imgFilesList.name})';
imgFileIndx = str2double(erase(imgFileNames, "img_" | ".bin"));
if ~strcmpi(imgFileNames{1}, 'img_00000.bin') || any(diff(imgFileIndx) ~= 1)
    error('Image binary files (img_xxxxx.bin) missing! Classification aborted.')
end

header = memmapfile([DataFolder imgFilesList(1).name], ...
    'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
nx = double(header.Data.header(2));
ny = double(header.Data.header(3));
frameFormat = {'uint64', 3, 'framej'; 'uint16', [double(nx), double(ny)], 'imgj'};
ImRes_XY = [nx, ny];
SizeImage = nx * ny * 2 + 3 * 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SubROI
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
            dat = memmapfile([DataFolder imgFilesList(1).name], ...
                'Offset', hWima * 4 + 5 * SizeImage, ...
                'Format', frameFormat, 'repeat', 1);
            dat = dat.Data.imgj;
            fig = figure;
            imagesc(dat);
            h = drawrectangle();
            wait(h);
            Pos = h.Position;
            close(fig);
        case 'Cancel'
            disp('User pressed cancel')
            Pos = [1 1 ImRes_XY(1) ImRes_XY(2)];
    end

    LimX = [round(Pos(1)) round(Pos(1)+Pos(3))];
    LimY = [round(Pos(2)) round(Pos(2)+Pos(4))];
    save([SaveFolder 'ROI.mat'], 'Pos');
    outFile{end+1} = 'ROI.mat'; %#ok<AGROW>
else
    LimX = [1 ImRes_XY(1)];
    LimY = [1 ImRes_XY(2)];
end

Rx = round((LimX(2) - LimX(1) + 1) / BinningSpatial);
Ry = round((LimY(2) - LimY(1) + 1) / BinningSpatial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How many colors and in which order?
fprintf('Sorting images per channels.\n');
if AcqInfoStream.MultiCam
    Tags = fieldnames(AcqInfoStream);
    idx = contains(Tags, 'Illumination');
    NbColors = sum(idx);
    Colors = struct('ID', {}, 'Color', {}, 'CamIdx', {}, 'FrameIdx', {}, 'Exposure', {});
    for indC = 1:NbColors
        Colors(indC).ID = indC;
        eval(['Colors(' int2str(indC) ').Color = AcqInfoStream.Illumination' int2str(indC) '.Color;']);
        eval(['Colors(' int2str(indC) ').CamIdx = AcqInfoStream.Illumination' int2str(indC) '.CamIdx;']);
        eval(['Colors(' int2str(indC) ').FrameIdx = AcqInfoStream.Illumination' int2str(indC) '.FrameIdx;']);
        Colors(indC).Exposure = localResolveExposure(AcqInfoStream, Colors(indC).Color);
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
    channelInfo = [channelInfo1, channelInfo2];
    fprintf('Camera #2 Done.\n');
else
    Tags = fieldnames(AcqInfoStream);
    idx = contains(Tags, 'Illumination');
    NbColors = sum(idx);
    Colors = struct('ID', {}, 'Color', {}, 'Exposure', {});
    for indC = 1:NbColors
        Colors(indC).ID = indC;
        eval(['Colors(' int2str(indC) ').Color = AcqInfoStream.Illumination' int2str(indC) '.Color;']);
        Colors(indC).Exposure = localResolveExposure(AcqInfoStream, Colors(indC).Color);
    end

    imgFilesList = dir([DataFolder 'img_*.bin']);
    channelInfo = ChannelsSort(imgFilesList, Colors);
    fprintf('Done.\n');
end

% Check if all color channels have the same length:
datLenList = [channelInfo.datLength];
if ~isscalar(unique(datLenList))
    disp('Fixing data length...')
    newLen = min(datLenList);
    idxLong = find(datLenList > newLen);
    for ind = idxLong
        datPath = fullfile(SaveFolder, channelInfo(ind).datFile);
        fid = fopen(datPath, 'r');
        dat = fread(fid, Inf, '*single');
        fclose(fid);
        dat = reshape(dat, channelInfo(ind).datSize(1), channelInfo(ind).datSize(2), []);
        dat = dat(:,:,1:newLen);
        fid = fopen(datPath, 'w');
        fwrite(fid, dat, 'single');
        fclose(fid);
        channelInfo(ind).datLength = newLen;
    end
    disp('Data length fixed.')
end

% Save shared acquisition metadata for the processed outputs.
AcqInfoStream.Width = Rx;
AcqInfoStream.Height = Ry;
AcqInfoStream.Length = channelInfo(1).datLength;
AcqInfoStream.FrameRateHz = AcqInfoStream.FrameRateHz / (numel(Colors) * BinningTemp);
AcqInfoStream.BinningSpatial = BinningSpatial;
AcqInfoStream.BinningTemp = BinningTemp;
save([SaveFolder 'AcqInfos.mat'], 'AcqInfoStream');
outFile = [outFile, {channelInfo.datFile}, {'AcqInfos.mat'}];

if nargout ~= 0
    varargout{1} = outFile;
end

    function channelInfoOut = ChannelsSort(fList, colors)
        % CHANNELSORT performs the classification of the raw data into the
        % existing colors and saves the data into separate .dat files.

        subNbColors = size(colors, 2);
        channelInfoOut = repmat(struct('datFile', '', 'datSize', [Ry, Rx], 'datLength', 0), 1, subNbColors);
        fid = zeros(1, subNbColors);

        for indC = 1:subNbColors
            dTag = localGetDatTag(colors(indC).Color);
            channelInfoOut(indC).datFile = dTag;
            channelInfoOut(indC).datSize = [Ry, Rx];
            channelInfoOut(indC).datLength = 0;
            if exist([SaveFolder dTag], 'file')
                delete([SaveFolder dTag]);
            end
            fid(indC) = fopen([SaveFolder dTag], 'w');
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
            overflow = mod(size(iData,3), subNbColors * BinningTemp);
            if overflow > 0
                oIm = iData(:,:, size(iData,3) - (overflow:-1:1) + 1);
            else
                oIm = [];
            end
            Images = iData(:,:,1:(size(iData,3) - overflow));
            clear iData hData overflow;

            if isempty(Images)
                break
            end
            Images = reshape(Images, ImRes_XY(2), ImRes_XY(1), subNbColors, []);

            for indC = 1:subNbColors
                Ims = squeeze(Images(:, :, indC, :));
                if any(sum(sum(Ims,1),2) == 0)
                    idxGood = find(sum(sum(Ims,1),2) > 1);
                    Ims = interp1(idxGood, single(reshape(Ims(:,:,idxGood),[], length(idxGood)))', 1:size(Ims,3), 'linear', 'extrap');
                    Ims = reshape(Ims', ImRes_XY(2), ImRes_XY(1), []);
                end

                if b_SubROI
                    Ims = Ims(round(LimY(1)):round(LimY(2)), round(LimX(1)):round(LimX(2)), :);
                end
                if BinningTemp > 1
                    Ims = imresize3(Ims, [size(Ims,1), size(Ims,2), size(Ims,3) / BinningTemp], 'linear');
                end
                if BinningSpatial > 1
                    Ims = imresize(Ims, 1 / BinningSpatial);
                end

                fwrite(fid(indC), single(Ims), 'single');
                channelInfoOut(indC).datLength = channelInfoOut(indC).datLength + size(Ims,3);
            end
            fprintf('\n');
        end

        for indC = 1:subNbColors
            fclose(fid(indC));
        end
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

function dTag = localGetDatTag(colorName)
%LOCALGETDATTAG Convert illumination color metadata to output .dat name.

if contains(colorName, {'red','green'}, 'IgnoreCase', true)
    dTag = [lower(colorName) '.dat'];
elseif contains(colorName, 'amber', 'IgnoreCase', true)
    dTag = 'yellow.dat';
elseif contains(colorName, 'fluo', 'IgnoreCase', true)
    waveTag = regexp(colorName, '[0-9]{3}', 'match');
    if ~isempty(waveTag)
        dTag = ['fluo_' waveTag{:} '.dat'];
    else
        dTag = 'fluo.dat';
    end
else
    dTag = 'speckle.dat';
end
end
