function outData = genRetinotopyMaps(data, SaveFolder, varargin)
%GENRETINOTOPYMAPS Generate retinotopy amplitude and phase maps.
%
%   outData = genRetinotopyMaps(data, SaveFolder)
%   outData = genRetinotopyMaps(data, SaveFolder, 'Name', Value, ...)
%
%   This function generates retinotopy azimuth and/or elevation maps from
%   an image time series and stimulation timing information stored in
%   events.mat.
%
%   Inputs:
%       data       - Numeric Y x X x T array, or raw .dat filename.
%       SaveFolder - Folder containing events.mat.
%
%   Name-Value parameters:
%       b_useAverageMovie - Logical scalar. If true, compute maps from the
%                           average trial movie. Default: false
%       Direction         - 'All', 'Azimuth_only', or 'Elevation_only'.
%                           Default: 'All'
%       ViewingDist_cm    - Viewing distance in cm. Default: 0
%       ScreenXsize_cm    - Screen width in cm. Default: 0
%       ScreenYsize_cm    - Screen height in cm. Default: 0
%       FrameRateHz       - Frame rate of DATA (Hz). PipelineManager injects
%                           it from the data; a .dat input's header provides
%                           it otherwise. Numeric input needs it explicitly;
%                           AcqInfos.mat is not used.
%
%   Output:
%       outData    - UMT struct with one or two entries:
%                    * AzimuthMap   : Y x X x F
%                    * ElevationMap : Y x X x F
%                    where F = {'Amplitude','Phase'}.
%
%   Notes:
%       - events.mat must contain: eventID, timestamps, state, eventNameList.
%       - Standard direction labels are '0', '90', '180', and '270'.
%       - If events.mat uses generic labels, the function can locally
%         standardize them according to Direction without modifying the
%         file on disk. A warning is raised when this assumption is used.
%       - With b_useAverageMovie = false, only ON-to-OFF stimulus epochs
%         are concatenated before the FFT. Inter-stimulus gaps are excluded.
%         This makes bin nKeptSweeps+1 represent the stimulus repetition in the
%         concatenated window, at the cost of discontinuities between epochs.
%         Previous releases transformed one contiguous first-onset-to-last-
%         offset span, so their maps are not numerically comparable.
%       - For raw .dat input, both FFT paths process spatial X slabs. The
%         final Y x X amplitude and phase maps remain resident in RAM, while
%         average-movie and baseline scratch storage scales with slab width.

% Default output for pipeline management.
default_Output = 'retinotopyMaps.umt';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) && ...
        strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = 'genRetinotopyMaps';
addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addParameter(p, 'b_useAverageMovie', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Direction', 'All', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ViewingDist_cm', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'ScreenXsize_cm', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'ScreenYsize_cm', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'FrameRateHz', []);
parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
opts = struct();
opts.b_useAverageMovie = p.Results.b_useAverageMovie;
opts.Direction = char(string(p.Results.Direction));
opts.ViewingDist_cm = double(p.Results.ViewingDist_cm);
opts.ScreenXsize_cm = double(p.Results.ScreenXsize_cm);
opts.ScreenYsize_cm = double(p.Results.ScreenYsize_cm);

% -------------------------------------------------------------------------
% Resolve input data and metadata
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty'}, mfilename, 'data');
    assert(ndims(data) == 3, ...
        'Umitoolbox:genRetinotopyMaps:WrongInput', ...
        'Numeric input must be a 3-D Y x X x T array.');
    dataIn = single(data);
    metaData = iResolveRepresentativeMeta(dataIn, ...
        resolveDataInfoValue('frameRateHz', p.Results.FrameRateHz, data, mfilename));
else
    assert(ischar(data) || (isstring(data) && isscalar(data)), ...
        'Umitoolbox:genRetinotopyMaps:WrongInput', ...
        ['Input must be a numeric Y x X x T array or a .dat filename. ' ...
         'UMT structs are not supported.']);
    dataFile = char(string(data));
    if ~isfile(dataFile)
        altPath = fullfile(SaveFolder, dataFile);
        if isfile(altPath)
            dataFile = altPath;
        else
            error('Umitoolbox:genRetinotopyMaps:MissingInput', ...
                'Input file "%s" was not found.', char(string(data)));
        end
    end

    [~,~,ext] = fileparts(dataFile);
    assert(strcmpi(ext, '.dat'), ...
        'Umitoolbox:genRetinotopyMaps:WrongInput', ...
        'Only numeric YXT input or raw .dat input are supported.');

    fileInfo = loadMetaData(dataFile);
    assertDatLayout(fileInfo, {{'Y','X','T'}}, 'genRetinotopyMaps');
    metaData = iInternalMetaFromInfo(fileInfo);
    if ~isempty(p.Results.FrameRateHz)
        % An explicit (or injected) rate wins over the header, with a
        % warning when they differ (resolveDataInfoValue).
        metaData.Freq = resolveDataInfoValue('frameRateHz', p.Results.FrameRateHz, ...
            dataFile, mfilename);
    end
    dataIn = dataFile;
end

assert(isequal(cellstr(string(metaData.dim_names(:).')), {'Y','X','T'}), ...
    'Umitoolbox:genRetinotopyMaps:WrongInput', ...
    'Input data must have dimensions Y, X, and T.');

assert(isfile(fullfile(SaveFolder, 'events.mat')), ...
    'Umitoolbox:genRetinotopyMaps:MissingInput', ...
    '"events.mat" file not found.');

evObj = EventsManager(SaveFolder);
% Ignored sweeps (selectedEvents false) keep their place in the sweep
% structure but are left out of the maps (.dat header Phase 8c).
selected = true(size(evObj.eventID));
if ~isempty(evObj.selectedEvents) && isequal(size(evObj.selectedEvents), size(evObj.eventID))
    selected = logical(evObj.selectedEvents);
end
evntInfo = struct( ...
    'eventID', evObj.eventID, ...
    'timestamps', evObj.timestamps, ...
    'state', evObj.state, ...
    'selected', selected, ...
    'eventNameList', {evObj.eventNameList});
evntInfo = iStandardizeEvents(evntInfo, opts.Direction);

% -------------------------------------------------------------------------
% Dispatch processing mode
% -------------------------------------------------------------------------
if ischar(dataIn) || (isstring(dataIn) && isscalar(dataIn))
    mapStruct = chunkedDatMode(char(string(dataIn)), metaData, evntInfo, opts);
else
    mapStruct = standardMode(dataIn, metaData, evntInfo, opts);
end

outData = iPackageOutputUMT(mapStruct);

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo( ...
            'genRetinotopyMaps', ...
            'Generate retinotopy azimuth and elevation maps from image time series.');

        info = PipelineManager.addInput(info, ...
            'data', ...
            'ImageTimeSeries', ...
            'Numeric YXT array or raw .dat filename.', ...
            'position', 1, ...
            'callType', 'positional', ...
            'isData', true, ...
            'supportsFile', true, ...
            'dataMode', 'either');

        info = PipelineManager.addInput(info, ...
            'SaveFolder', ...
            'SaveFolder', ...
            'Folder containing events.mat and metadata sources.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput(info, ...
            'b_useAverageMovie', ...
            'parameter', ...
            'If true, compute maps from the average trial movie.', ...
            'kind', 'parameter', ...
            'position', 3, ...
            'callType', 'namevalue', ...
            'default', false, ...
            'allowed', [true false], ...
            'dataType', 'logical');

        info = PipelineManager.addInput(info, ...
            'Direction', ...
            'parameter', ...
            'Expected stimulus directions.', ...
            'kind', 'parameter', ...
            'position', 4, ...
            'callType', 'namevalue', ...
            'default', 'All', ...
            'allowed', {'All','Azimuth_only','Elevation_only'}, ...
            'dataType', 'char');

        info = PipelineManager.addInput(info, ...
            'ViewingDist_cm', ...
            'parameter', ...
            'Viewing distance in cm.', ...
            'kind', 'parameter', ...
            'position', 5, ...
            'callType', 'namevalue', ...
            'default', 0, ...
            'allowed', [0 Inf], ...
            'dataType', 'numeric');

        info = PipelineManager.addInput(info, ...
            'ScreenXsize_cm', ...
            'parameter', ...
            'Screen width in cm.', ...
            'kind', 'parameter', ...
            'position', 6, ...
            'callType', 'namevalue', ...
            'default', 0, ...
            'allowed', [0 Inf], ...
            'dataType', 'numeric');

        info = PipelineManager.addInput(info, ...
            'ScreenYsize_cm', ...
            'parameter', ...
            'Screen height in cm.', ...
            'kind', 'parameter', ...
            'position', 7, ...
            'callType', 'namevalue', ...
            'default', 0, ...
            'allowed', [0 Inf], ...
            'dataType', 'numeric');

        info = PipelineManager.addInput(info, ...
            'FrameRateHz', ...
            'sourceInfo', ...
            'Frame rate of the input data (Hz), injected from the data.', ...
            'kind', 'sourceInfo', ...
            'sourceField', 'frameRateHz', ...
            'required', false);

        info = PipelineManager.addOutput(info, ...
            'outData', ...
            'ProcessedData', ...
            'data', ...
            'Retinotopy output UMT with AzimuthMap and/or ElevationMap entries.', ...
            default_Output, ...
            1, ...
            'isData', true);
    end
end

%% ==================== STANDARD MODE ====================
function mapStruct = standardMode(data, metaData, evntInfo, opts)
ampMaps = cell(size(evntInfo.eventNameList));
phiMaps = ampMaps;

framestamps = round(evntInfo.timestamps * metaData.Freq);
iValidateEventTransitions(evntInfo, framestamps, size(data, 3), opts.b_useAverageMovie);
w = waitbar(0,'Calculating FFT ...','Name','genRetinotopyMaps');

for ind = 1:numel(evntInfo.eventNameList)
    waitbar((ind-1)/numel(evntInfo.eventNameList), w, ...
        ['Calculating FFT for direction ' strrep(evntInfo.eventNameList{ind}, '_', ' ')]);
    [indxOnAll, indxOffAll, indxOn, indxOff] = iDirectionSweeps(evntInfo, ind);
    if isempty(indxOn)
        % Every sweep of this direction is ignored: NaN maps (Phase 8c).
        ampMaps{ind} = nan(metaData.datSize, 'single');
        phiMaps{ind} = nan(metaData.datSize, 'single');
        waitbar(ind/numel(evntInfo.eventNameList),w);
        continue
    end
    freqFFT = iFrequencyIndex(opts, numel(indxOn));

    if opts.b_useAverageMovie
        % Window lengths come from every sweep; only kept sweeps are averaged.
        bsln_len = round(mean(framestamps(indxOnAll(2:end)) - framestamps(indxOffAll(1:end-1))));
        trial_len = round(mean(framestamps(indxOffAll) - framestamps(indxOnAll)));
        assert(~isempty(bsln_len) && bsln_len > 0, ...
            'Umitoolbox:genRetinotopyMaps:InvalidBaseline', ...
            ['Could not compute a valid baseline period. Check whether the ' ...
             'recording contains inter-stimulus time or disable b_useAverageMovie.']);
        assert(trial_len >= freqFFT, ...
            'Umitoolbox:genRetinotopyMaps:EpochsTooShort', ...
            ['Direction %s has %d-frame stimulus epochs, fewer than the %d needed ' ...
             'to resolve the stimulus frequency. Check the frame rate and the ' ...
             'epoch timestamps in events.mat.'], ...
            evntInfo.eventNameList{ind}, trial_len, freqFFT);

        avg_mov = zeros([metaData.datSize, trial_len + bsln_len],'single');
        for ii = 1:length(indxOn)
            tStart = framestamps(indxOn(ii)) - bsln_len;
            tEnd   = framestamps(indxOn(ii)) + trial_len - 1;
            assert(tStart >= 1 && tEnd <= size(data, 3), ...
                'Umitoolbox:genRetinotopyMaps:InvalidBaseline', ...
                ['Could not fit the average-movie baseline and trial windows ' ...
                 'within the recording for trial %d of direction %d.'], ii, ind);

            DeltaR = data(:,:,tStart:tEnd) - ...
                median(data(:,:,tStart:framestamps(indxOn(ii))-1), 3,'omitnan');
            avg_mov = avg_mov + DeltaR;
        end
        avg_mov = avg_mov / length(indxOn);
        fDat = fft(avg_mov(:,:,bsln_len+1:end),[],3);
    else
        frOn = [];
        for ii = 1:length(indxOn)
            frOn = [frOn, framestamps(indxOn(ii)):framestamps(indxOff(ii))-1]; %#ok<AGROW>
        end
        assert(numel(frOn) >= freqFFT, ...
            'Umitoolbox:genRetinotopyMaps:EpochsTooShort', ...
            ['Direction %s has %d stimulus frame(s) in its kept sweeps, fewer than ' ...
             'the %d needed to resolve the stimulus frequency. Check the frame rate ' ...
             'and the epoch timestamps in events.mat.'], ...
            evntInfo.eventNameList{ind}, numel(frOn), freqFFT);
        fDat = fft(data(:,:,frOn),[],3);
    end

    ampMaps{ind} = (abs(fDat(:,:,freqFFT)) * 2) / size(fDat,3);
    phiMaps{ind} = mod(-angle(fDat(:,:,freqFFT)),2*pi);
    waitbar(ind/numel(evntInfo.eventNameList),w);
end
close(w);

mapStruct = buildMaps(ampMaps, phiMaps, metaData, opts, evntInfo);
end

%% ==================== CHUNKED RAW-DAT MODE ====================
function mapStruct = chunkedDatMode(datFile, metaData, evntInfo, opts)
nY = metaData.datSize(1);
nX = metaData.datSize(2);
nT = metaData.datLength;

ampMaps = cell(size(evntInfo.eventNameList));
phiMaps = ampMaps;
framestamps = round(evntInfo.timestamps * metaData.Freq);
iValidateEventTransitions(evntInfo, framestamps, nT, opts.b_useAverageMovie);


w = waitbar(0,'Calculating FFT (Low RAM usage) ...','Name','genRetinotopyMaps');
slabIn = spatialSlabIO('open', datFile, 'Info', metaData);
cIn = onCleanup(@() spatialSlabIO('close', slabIn));

for ind = 1:numel(evntInfo.eventNameList)
    waitbar((ind-1)/numel(evntInfo.eventNameList), w, ...
        ['Calculating FFT for direction ' strrep(evntInfo.eventNameList{ind}, '_', ' ')]);
    [indxOnAll, indxOffAll, indxOn, indxOff] = iDirectionSweeps(evntInfo, ind);
    if isempty(indxOn)
        % Every sweep of this direction is ignored: NaN maps (Phase 8c).
        ampMaps{ind} = nan(nY,nX,'single');
        phiMaps{ind} = nan(nY,nX,'single');
        waitbar(ind/numel(evntInfo.eventNameList),w);
        continue
    end
    freqFFT = iFrequencyIndex(opts, numel(indxOn));

    ampMaps{ind} = zeros(nY,nX,'single');
    phiMaps{ind} = zeros(nY,nX,'single');

    if opts.b_useAverageMovie
        % Window lengths come from every sweep; only kept sweeps are averaged.
        bsln_len = round(mean(framestamps(indxOnAll(2:end)) - framestamps(indxOffAll(1:end-1))));
        trial_len = round(mean(framestamps(indxOffAll) - framestamps(indxOnAll)));
        assert(~isempty(bsln_len) && bsln_len > 0, ...
            'Umitoolbox:genRetinotopyMaps:InvalidBaseline', ...
            ['Could not compute a valid baseline period. Check whether the ' ...
             'recording contains inter-stimulus time or disable b_useAverageMovie.']);
        assert(trial_len >= freqFFT, ...
            'Umitoolbox:genRetinotopyMaps:EpochsTooShort', ...
            ['Direction %s has %d-frame stimulus epochs, fewer than the %d needed ' ...
             'to resolve the stimulus frequency. Check the frame rate and the ' ...
             'epoch timestamps in events.mat.'], ...
            evntInfo.eventNameList{ind}, trial_len, freqFFT);

        total_len = trial_len + bsln_len;

        % Estimate all slab-sized live arrays: average accumulator,
        % baseline frames, complex FFT output, and small frame/map scratch.
        bytesPerX = double(nY) * 4 * ...
            (double(total_len) + double(bsln_len) + ...
             2 * double(trial_len) + 4);
        nChunks = calculateMaxChunkSize(bytesPerX * double(nX), 1, .1);
        chunkX = max(1, ceil(nX / nChunks));
        nChunks = ceil(nX / chunkX);

        ampMap = ampMaps{ind};
        phiMap = phiMaps{ind};
        for c = 1:nChunks
            xStart = (c-1) * chunkX + 1;
            xEnd = min(xStart + chunkX - 1, nX);
            xIdx = xStart:xEnd;

            avgSlab = zeros(nY, numel(xIdx), total_len, 'single');
            for ii = 1:length(indxOn)
                tStart = framestamps(indxOn(ii)) - bsln_len;
                tEnd = framestamps(indxOn(ii)) + trial_len - 1;
                assert(tStart >= 1 && tEnd <= nT, ...
                    'Umitoolbox:genRetinotopyMaps:InvalidBaseline', ...
                    ['Could not fit the average-movie baseline and trial windows ' ...
                     'within the recording for trial %d of direction %d.'], ii, ind);

                baselineData = iReadFramesXSlab(slabIn, tStart:(tStart + bsln_len - 1), xIdx);
                baseline = median(baselineData, 3, 'omitnan');
                avgSlab(:,:,1:bsln_len) = avgSlab(:,:,1:bsln_len) + ...
                    (baselineData - baseline);

                trialData = iReadFramesXSlab(slabIn, (tStart + bsln_len):tEnd, xIdx);
                avgSlab(:,:,(bsln_len + 1):total_len) = ...
                    avgSlab(:,:,(bsln_len + 1):total_len) + (trialData - baseline);
            end

            avgSlab = avgSlab / length(indxOn);
            fSlab = fft(avgSlab(:,:,bsln_len+1:end), [], 3);
            ampMap(:,xIdx) = (abs(fSlab(:,:,freqFFT)) * 2) / size(fSlab,3);
            phiMap(:,xIdx) = mod(-angle(fSlab(:,:,freqFFT)), 2*pi);
        end
        ampMaps{ind} = ampMap;
        phiMaps{ind} = phiMap;
        waitbar(ind/numel(evntInfo.eventNameList),w);

    else
        frOn = [];
        for ii = 1:length(indxOn)
            frOn = [frOn, framestamps(indxOn(ii)):framestamps(indxOff(ii))-1]; %#ok<AGROW>
        end
        assert(numel(frOn) >= freqFFT, ...
            'Umitoolbox:genRetinotopyMaps:EpochsTooShort', ...
            ['Direction %s has %d stimulus frame(s) in its kept sweeps, fewer than ' ...
             'the %d needed to resolve the stimulus frequency. Check the frame rate ' ...
             'and the epoch timestamps in events.mat.'], ...
            evntInfo.eventNameList{ind}, numel(frOn), freqFFT);

        ampMap = ampMaps{ind};
        phiMap = phiMaps{ind};

        % Only this direction's stimulus frames are read (sorted and unique,
        % as spatialSlabIO requires); frMap restores the sweep order and any
        % repeated frame of frOn. The slab holds the frames, their complex
        % FFT (two real copies), and a reordered copy when frMap is needed.
        [frUnique, ~, frMap] = unique(frOn);
        frMap = frMap(:).';
        bReorder = ~isequal(frMap, 1:numel(frOn));
        nChunks = calculateMaxChunkSize(nY * nX * numel(frUnique) * 4, 4, .1);
        chunkX  = ceil(nX / nChunks);
        nChunks = ceil(nX / chunkX);

        for c = 1:nChunks
            fprintf('Processing spatial slab [%i/%i] for direction %s\n', c, nChunks, evntInfo.eventNameList{ind})
            xStart = (c-1)*chunkX + 1;
            xEnd   = min(xStart + chunkX - 1, nX);
            xIdx   = xStart:xEnd;

            slabData = iReadFramesXSlab(slabIn, frUnique, xIdx);
            if bReorder
                slabData = slabData(:,:,frMap);
            end

            fSlab = fft(slabData,[],3);
            ampSlab = (abs(fSlab(:,:,freqFFT))*2)/size(fSlab,3);
            phiSlab = mod(-angle(fSlab(:,:,freqFFT)),2*pi);
            ampMap(:,xIdx) = ampSlab;
            phiMap(:,xIdx) = phiSlab;
        end
        ampMaps{ind} = ampMap;
        phiMaps{ind} = phiMap;
        waitbar(ind/numel(evntInfo.eventNameList),w);
    end
end
close(w);

mapStruct = buildMaps(ampMaps, phiMaps, metaData, opts, evntInfo);
end

function frames = iReadFramesXSlab(slabIn, frameIdx, xIdx)
%IREADFRAMESXSLAB Read X columns of consecutive frames through spatialSlabIO.

frames = single(spatialSlabIO('read', slabIn, xIdx, frameIdx));
frames = reshape(frames, size(frames, 1), numel(xIdx), numel(frameIdx));
end

function metaData = iInternalMetaFromInfo(metaData)
%IINTERNALMETAFROMINFO Add this function's internal fields to loadMetaData output.
%
% dim_names, datSize = [Y X], datLength, and Freq are the fields used
% throughout this function (the same model iResolveRepresentativeMeta
% builds for numeric input). They are derived from the .dat Info schema
% only; the schema fields are kept for spatialSlabIO('open', ..., 'Info').

metaData.dim_names = metaData.dimNames;
metaData.datSize = [datAxisSize(metaData, 'Y'), datAxisSize(metaData, 'X')];
metaData.datLength = datAxisSize(metaData, 'T');
metaData.Freq = double(metaData.frameRateHz);
end

%% ==================== BUILD MAPS ====================
function mapStruct = buildMaps(ampMaps, phiMaps, metaData, opts, evntInfo)
mapStruct = struct();
[idxAz, indxAz] = ismember({'0','180'}, evntInfo.eventNameList);
[idxEl, indxEl] = ismember({'90','270'}, evntInfo.eventNameList);

if all(idxAz)
    azimMap = zeros([metaData.datSize 2], 'single');
    azimMap(:,:,1) = mean(cat(3, ampMaps{indxAz}), 3);
    azimMap(:,:,2) = pi + ((phiMaps{indxAz(1)} - phiMaps{indxAz(2)})/2);
    mapStruct.AzimuthMap = azimMap;
end

if all(idxEl)
    elevMap = zeros([metaData.datSize 2], 'single');
    elevMap(:,:,1) = mean(cat(3, ampMaps{indxEl}), 3);
    elevMap(:,:,2) = pi + ((phiMaps{indxEl(1)} - phiMaps{indxEl(2)})/2);
    mapStruct.ElevationMap = elevMap;
end

if all([opts.ViewingDist_cm, opts.ScreenXsize_cm, opts.ScreenYsize_cm] > 0)
    if isfield(mapStruct, 'AzimuthMap')
        va_az = atand(opts.ScreenXsize_cm/(2*opts.ViewingDist_cm));
        va_az = [-va_az va_az];
        mapStruct.AzimuthMap(:,:,2) = rescale(mapStruct.AzimuthMap(:,:,2), va_az(1), va_az(2));
    end
    if isfield(mapStruct, 'ElevationMap')
        va_el = atand(opts.ScreenYsize_cm/(2*opts.ViewingDist_cm));
        va_el = [-va_el va_el];
        mapStruct.ElevationMap(:,:,2) = rescale(mapStruct.ElevationMap(:,:,2), va_el(1), va_el(2));
    end
elseif any([opts.ViewingDist_cm, opts.ScreenXsize_cm, opts.ScreenYsize_cm] > 0)
    warning('Umitoolbox:genRetinotopyMaps:VisualAngleIncomplete', ...
        'Cannot calculate visual angle. Phase values will remain in radians.');
end
end

%% ==================== EVENT STANDARDIZATION ====================
function [indxOnAll, indxOffAll, indxOn, indxOff] = iDirectionSweeps(evntInfo, ind)
%IDIRECTIONSWEEPS ON/OFF transitions of direction IND: all, and kept ones.
%   Ignored sweeps (evntInfo.selected false) keep their place for timing
%   estimates but are left out of the maps (.dat header Phase 8c).
indxOnAll = find(evntInfo.eventID == ind & evntInfo.state == 1);
indxOffAll = find(evntInfo.eventID == ind & evntInfo.state == 0);
if isfield(evntInfo, 'selected')
    indxOn = indxOnAll(evntInfo.selected(indxOnAll));
    indxOff = indxOffAll(evntInfo.selected(indxOffAll));
else
    indxOn = indxOnAll;
    indxOff = indxOffAll;
end
end

function freqFFT = iFrequencyIndex(opts, nKeptSweeps)
%IFREQUENCYINDEX FFT bin of the stimulus frequency.
%   Average movie: one sweep per window (bin 2). Concatenated sweeps: one
%   cycle per kept sweep (bin nKeptSweeps + 1).
if opts.b_useAverageMovie
    freqFFT = 2;
else
    freqFFT = nKeptSweeps + 1;
end
end

function evntInfo = iStandardizeEvents(evntInfo, directionMode)
requiredFields = {'eventID','timestamps','state','eventNameList'};
assert(all(isfield(evntInfo, requiredFields)), ...
    'Umitoolbox:genRetinotopyMaps:MissingInput', ...
    'events.mat must contain eventID, timestamps, state, and eventNameList.');

switch lower(string(directionMode))
    case "all"
        nDir = 4;
        evNames = {'0', '90', '180', '270'};
    case "azimuth_only"
        nDir = 2;
        evNames = {'0', '180'};
    case "elevation_only"
        nDir = 2;
        evNames = {'90', '270'};
    otherwise
        error('Umitoolbox:genRetinotopyMaps:InvalidDirection', ...
            'Unknown Direction option: %s.', string(directionMode));
end

nOnsets = sum(evntInfo.state == 1);
nSweeps = nOnsets / nDir;

assert(abs(nSweeps - round(nSweeps)) < eps(max(1, nSweeps)), ...
    'Umitoolbox:genRetinotopyMaps:InvalidSweeps', ...
    ['The number of active stimulus epochs is not compatible with Direction = "%s". ' ...
     'The inferred number of sweeps is not an integer. Please check events.mat.'], ...
    string(directionMode));

evntInfo.nSweeps = round(nSweeps);

hasStandardLabels = any(ismember(evntInfo.eventNameList, {'0','90','180','270'}));
if ~hasStandardLabels
    warning('Umitoolbox:genRetinotopyMaps:GenericDirectionLabels', ...
        ['events.mat does not contain standard direction labels (0, 90, 180, 270). ' ...
         'Assuming direction order "%s" with %d sweep(s) and standardizing labels locally.'], ...
        strjoin(evNames, ', '), evntInfo.nSweeps);

    onIdx = find(evntInfo.state == 1);
    offIdx = find(evntInfo.state == 0);
    assert(numel(onIdx) == numel(offIdx), ...
        'Umitoolbox:genRetinotopyMaps:InvalidEvents', ...
        'events.mat must contain matching numbers of onset and offset timestamps.');

    pairedIDs = repmat((1:nDir)', evntInfo.nSweeps, 1);
    assert(numel(pairedIDs) == numel(onIdx), ...
        'Umitoolbox:genRetinotopyMaps:InvalidSweeps', ...
        ['The number of active stimulus epochs is not compatible with Direction = "%s". ' ...
         'The inferred number of sweeps is not an integer. Please check events.mat.'], ...
        string(directionMode));

    newEventID = zeros(size(evntInfo.eventID));
    newEventID(onIdx) = pairedIDs;
    newEventID(offIdx) = pairedIDs;

    evntInfo.eventNameList = evNames;
    evntInfo.eventID = newEventID;
end

assert(isequal(sort(string(evntInfo.eventNameList(:))), sort(string(evNames(:)))), ...
    'Umitoolbox:genRetinotopyMaps:InvalidDirectionNames', ...
    'Invalid direction names in events.mat. Expected labels: %s.', strjoin(evNames, ', '));
end

function iValidateEventTransitions(evntInfo, framestamps, nFrames, bUseAverageMovie)
%IVALIDATEEVENTTRANSITIONS Validate onset/offset pairing and frame bounds.

for ind = 1:numel(evntInfo.eventNameList)
    indxOn = find(evntInfo.eventID == ind & evntInfo.state == 1);
    indxOff = find(evntInfo.eventID == ind & evntInfo.state == 0);

    assert(~isempty(indxOn) && numel(indxOn) == numel(indxOff), ...
        'Umitoolbox:genRetinotopyMaps:InvalidEvents', ...
        ['Each retinotopy direction must contain the same non-zero number ' ...
         'of onset and offset timestamps. Direction %d is invalid.'], ind);

    onsetFrames = framestamps(indxOn);
    offsetFrames = framestamps(indxOff);
    assert(all(isfinite([onsetFrames(:); offsetFrames(:)])) && ...
        all(onsetFrames >= 1) && all(offsetFrames <= nFrames) && ...
        all(offsetFrames > onsetFrames), ...
        'Umitoolbox:genRetinotopyMaps:InvalidEvents', ...
        ['Retinotopy event transitions for direction %d must fall within ' ...
         'the recording and each offset must follow its onset.'], ind);

    if bUseAverageMovie
        assert(numel(indxOn) >= 2, ...
            'Umitoolbox:genRetinotopyMaps:InvalidBaseline', ...
            ['Average-movie mode requires at least two trials per direction ' ...
             'to estimate the inter-stimulus baseline. Direction %d is invalid.'], ind);

        baselineLength = round(mean(onsetFrames(2:end) - offsetFrames(1:end-1)));
        trialLength = round(mean(offsetFrames - onsetFrames));
        assert(isfinite(baselineLength) && baselineLength > 0 && ...
            isfinite(trialLength) && trialLength > 0 && ...
            onsetFrames(1) - baselineLength >= 1 && ...
            onsetFrames(end) + trialLength - 1 <= nFrames, ...
            'Umitoolbox:genRetinotopyMaps:InvalidBaseline', ...
            ['Could not fit the average-movie baseline and trial windows ' ...
             'within the recording for direction %d.'], ind);
    end
end
end

%% ==================== REPRESENTATIVE METADATA ====================
function metaData = iResolveRepresentativeMeta(data, frameRateHz)
%IRESOLVEREPRESENTATIVEMETA Build metadata for a raw numeric YXT array.
%
% datSize/datLength/dim_names come from the array being processed; the
% frame rate is the data's own (explicit or injected FrameRateHz; .dat
% header Phase 7a: AcqInfos.mat is not used).

metaData = struct();
metaData.dim_names = {'Y','X','T'};
metaData.datSize = [size(data,1), size(data,2)];
metaData.datLength = size(data,3);
metaData.Freq = double(frameRateHz);
end

%% ==================== PACKAGE OUTPUT ====================
function outUMT = iPackageOutputUMT(mapStruct)
assert(~isempty(fieldnames(mapStruct)), ...
    'Umitoolbox:genRetinotopyMaps:NoMapsGenerated', ...
    'No retinotopy maps could be generated from the provided directions.');

labels = struct();
labels.F = {'Amplitude', 'Phase'};

entryNames = fieldnames(mapStruct);
outUMT = [];
for iEntry = 1:numel(entryNames)
    if iEntry == 1
        outUMT = genUMTStruct( ...
            mapStruct.(entryNames{iEntry}), ...
            'kind', 'image', ...
            'entryName', entryNames{iEntry}, ...
            'dimNames', {'Y','X','F'}, ...
            'labels', labels);
    else
        outUMT = genUMTStruct( ...
            outUMT, ...
            'value', mapStruct.(entryNames{iEntry}), ...
            'entryName', entryNames{iEntry}, ...
            'dimNames', {'Y','X','F'});
    end
end
end
