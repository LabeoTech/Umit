classdef TestStreamingReadersHeaderedInput < matlab.unittest.TestCase
    %TESTSTREAMINGREADERSHEADEREDINPUT Streaming Analysis readers give identical
    %results on headerless and headered input.
    %
    %   .dat header Phases 3c-1 and 4b. Each function runs twice: once on a
    %   synthetic headerless folder (legacy sidecars since Phase 5a) and once
    %   on a copy in which every .dat
    %   carries a header (same values, same AcqInfos.mat and events.mat).
    %   Outputs are compared after loading file outputs with loadData and
    %   replacing folder paths in text. Synthetic folders follow
    %   TestAnalysisPipelineManagerIntegration and, for HemoCompute,
    %   TestHemoCompute. The Phase 4b cases cover the IOIAnalysis cores'
    %   input reads; the event-split NormalisationFiltering path has no
    %   headerless Y,X,T,E form and is covered by
    %   TestNormalisationFiltering.testEventSplitFileModeReturnsArray.
    %   Phases 4c-2a and 4c-2b: outputIsHeaderedWithInputInfo checks that
    %   the streamed Analysis writers and the IOIAnalysis cores write
    %   headered outputs that inherit the input's frame rate and exposure
    %   and carry the output name. For HemoCompute the input is the green
    %   channel, which is on the target (lowest-frequency) timeline.

    properties
        TempRoot
    end

    properties (TestParameter)
        caseSpec = struct( ...
            'GSR', struct('kind', 'image', ...
                'call', @(f) GSR('input.dat', f, 'UseMask', false)), ...
            'normalizeZScore', struct('kind', 'image', ...
                'call', @(f) normalizeZScore('input.dat', f)), ...
            'normalizeBSLN', struct('kind', 'image', ...
                'call', @(f) normalizeBSLN('input.dat', f)), ...
            'apply_detrend', struct('kind', 'image', ...
                'call', @(f) apply_detrend('input.dat', f)), ...
            'apply_aggregate_function', struct('kind', 'image', ...
                'call', @(f) apply_aggregate_function('input.dat', f)), ...
            'spatialGaussFilt', struct('kind', 'image', ...
                'call', @(f) spatialGaussFilt('input.dat', f)), ...
            'genRetinotopyMaps', struct('kind', 'retinotopy', ...
                'call', @(f) genRetinotopyMaps('input.dat', f)), ...
            'run_HemoCorrection', struct('kind', 'hemo', ...
                'call', @(f) run_HemoCorrection('yellow.dat', f, 'Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', false, 'Amber', false)), ...
            'run_BloodFlow_spatial', struct('kind', 'speckle', ...
                'call', @(f) run_BloodFlow(f, 'speckle.dat', 'sType', 'Spatial', 'ExposureMsec', 2)), ...
            'run_BloodFlow_temporal', struct('kind', 'speckle', ...
                'call', @(f) run_BloodFlow(f, 'speckle.dat', 'sType', 'Temporal', 'ExposureMsec', 2)), ...
            'HemoCompute_RAMsafe', struct('kind', 'hemocompute', ...
                'call', @(f) iHemoComputeOutputs(f)), ...
            'HemoCorrection_LinearRegression', struct('kind', 'hemo', ...
                'call', @(f) run_HemoCorrection('yellow.dat', f, 'Algorithm', 'LinearRegression', ...
                'Red', true, 'Green', false, 'Amber', false)), ...
            'NormalisationFiltering_YXT', struct('kind', 'image', ...
                'call', @(f) normalizeLPF('input.dat', f)), ...
            'Ana_Speckle_RAMsafe', struct('kind', 'speckle', ...
                'call', @(f) run_Ana_Speckle(f, 'speckle.dat')), ...
            'SpeckleMapping_RAMsafe_spatial', struct('kind', 'speckle', ...
                'call', @(f) run_SpeckleMapping(f, 'speckle.dat', 'sType', 'Spatial')), ...
            'SpeckleMapping_RAMsafe_temporal', struct('kind', 'speckle', ...
                'call', @(f) run_SpeckleMapping(f, 'speckle.dat', 'sType', 'Temporal')))

        writerCase = struct( ...
            'GSR', struct('kind', 'image', 'input', 'input.dat', ...
                'call', @(f) GSR('input.dat', f, 'UseMask', false)), ...
            'normalizeZScore', struct('kind', 'image', 'input', 'input.dat', ...
                'call', @(f) normalizeZScore('input.dat', f)), ...
            'spatialGaussFilt', struct('kind', 'image', 'input', 'input.dat', ...
                'call', @(f) spatialGaussFilt('input.dat', f)), ...
            'apply_detrend', struct('kind', 'image', 'input', 'input.dat', ...
                'call', @(f) apply_detrend('input.dat', f)), ...
            'run_HemoCorrection_ratiometric', struct('kind', 'hemo', 'input', 'yellow.dat', ...
                'call', @(f) run_HemoCorrection('yellow.dat', f, 'Algorithm', 'Ratiometric', ...
                'Red', true, 'Green', false, 'Amber', false)), ...
            'run_BloodFlow_spatial_normalized', struct('kind', 'speckle', 'input', 'speckle.dat', ...
                'call', @(f) run_BloodFlow(f, 'speckle.dat', 'sType', 'Spatial', ...
                'ExposureMsec', 2, 'bNormalize', true)), ...
            'run_BloodFlow_temporal', struct('kind', 'speckle', 'input', 'speckle.dat', ...
                'call', @(f) run_BloodFlow(f, 'speckle.dat', 'sType', 'Temporal', 'ExposureMsec', 2)), ...
            'Ana_Speckle_RAMsafe', struct('kind', 'speckle', 'input', 'speckle.dat', ...
                'call', @(f) run_Ana_Speckle(f, 'speckle.dat')), ...
            'Ana_Speckle_standard_saved', struct('kind', 'speckle', 'input', 'speckle.dat', ...
                'call', @(f) iAnaSpeckleStandardSaved(f)), ...
            'HemoCompute_RAMsafe', struct('kind', 'hemocompute', 'input', 'green.dat', ...
                'call', @(f) iHemoComputeHbOPath(f)), ...
            'HemoCorrection_LinearRegression', struct('kind', 'hemo', 'input', 'yellow.dat', ...
                'channelName', 'fluoHemoCorr', ...
                'call', @(f) run_HemoCorrection('yellow.dat', f, 'Algorithm', 'LinearRegression', ...
                'Red', true, 'Green', false, 'Amber', false)), ...
            'NormalisationFiltering_YXT', struct('kind', 'image', 'input', 'input.dat', ...
                'call', @(f) normalizeLPF('input.dat', f)))
    end

    methods (TestMethodSetup)
        function createTempRoot(testCase)
            testCase.TempRoot = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end
    end

    methods (Test)
        function headeredInputGivesIdenticalOutput(testCase, caseSpec)
            plainFolder = iCreateFolder(fullfile(testCase.TempRoot, 'plain'), caseSpec.kind);
            if strcmp(caseSpec.kind, 'hemocompute')
                rigStore = iCreateConfiguredRig();
                testCase.addTeardown(@() iRemoveFolder(rigStore.RigRoot));
                iBindDatasetToRig(plainFolder, rigStore);
            end
            headerFolder = iHeaderedCopy(plainFolder, fullfile(testCase.TempRoot, 'header'));

            testCase.assertFalse(iAnyHeadered(plainFolder), 'plain inputs must be headerless');
            testCase.assertTrue(iAllHeadered(headerFolder), 'copied inputs must be headered');

            plainOut = iOutputValue(caseSpec.call(plainFolder), plainFolder);
            headerOut = iOutputValue(caseSpec.call(headerFolder), headerFolder);

            testCase.verifyEqual(iNormalize(headerOut, headerFolder), ...
                iNormalize(plainOut, plainFolder));
        end

        function outputIsHeaderedWithInputInfo(testCase, writerCase)
            plainFolder = iCreateFolder(fullfile(testCase.TempRoot, 'plain'), writerCase.kind);
            if strcmp(writerCase.kind, 'hemocompute')
                rigStore = iCreateConfiguredRig();
                testCase.addTeardown(@() iRemoveFolder(rigStore.RigRoot));
                iBindDatasetToRig(plainFolder, rigStore);
            end
            headerFolder = iHeaderedCopy(plainFolder, fullfile(testCase.TempRoot, 'header'));

            folders = {plainFolder, headerFolder};
            for k = 1:numel(folders)
                folder = folders{k};
                inInfo = loadMetaData(fullfile(folder, writerCase.input));
                outFile = iOutputPath(writerCase.call(folder), folder);

                testCase.assertTrue(isDatWithHeader(outFile), ...
                    sprintf('%s must be headered', outFile));
                hdr = readDatHeader(outFile);
                [~, outBase] = fileparts(outFile);
                testCase.verifyEqual(hdr.frameRateHz, double(single(inInfo.frameRateHz)));
                testCase.verifyEqual(hdr.exposureMsec, double(single(inInfo.exposureMsec)));
                % run_HemoCorrection moves HemoCorrection's fluoHemoCorr.dat onto
                % its own output name; a move keeps the core's channelName.
                expectedName = outBase;
                if isfield(writerCase, 'channelName')
                    expectedName = writerCase.channelName;
                end
                testCase.verifyEqual(hdr.channelName, expectedName);
                testCase.verifyTrue(hdr.writeComplete);
                testCase.verifyEqual(hdr.dimSizes, inInfo.dimSizes);
            end
        end
    end
end

% =========================================================================
% Synthetic folders
% =========================================================================
function folder = iCreateFolder(folder, kind)
mkdir(folder);
if strcmp(kind, 'hemocompute')
    iCreateHemoComputeDataset(folder);
    return
end
Ny = 8;
Nx = 7;
frameRate = 10;
rng(17);

switch kind
    case 'image'
        Nt = 40;
        files = {'input.dat'};
        data = {randn(Ny, Nx, Nt, 'single')};
    case 'retinotopy'
        Nt = 140;
        files = {'input.dat'};
        data = {randn(Ny, Nx, Nt, 'single')};
    case 'hemo'
        Nt = 40;
        files = {'yellow.dat', 'red.dat'};
        data = {1 + 0.05 .* randn(Ny, Nx, Nt, 'single'), 1 + 0.05 .* randn(Ny, Nx, Nt, 'single')};
    case 'speckle'
        Nt = 40;
        files = {'speckle.dat'};
        data = {1 + 0.05 .* randn(Ny, Nx, Nt, 'single')};
end

for k = 1:numel(files)
    iWriteHeaderless(fullfile(folder, files{k}), data{k}, frameRate, 5);
end

AcqInfoStream = struct('Height', Ny, 'Width', Nx, 'Length', Nt, ...
    'FrameRateHz', frameRate, 'Freq', frameRate, 'Datatype', 'single', ...
    'ExposureMsec', 5, 'Camera_Model', 'SyntheticHeaderTest', 'MultiCam', false, ...
    'AISampleRate', 1000, 'AINChannels', 1, 'AICh1', 'CameraTrig');
save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
ensurePMReadyAcqInfos(folder, files);

switch kind
    case 'retinotopy'
        iWriteRetinotopyEvents(folder);
end
end

function headerFolder = iHeaderedCopy(plainFolder, headerFolder)
% Copy every file; rewrite each .dat with a header in front of the same
% data. The header takes each file's own axes, class, frame rate, and
% exposure from loadMetaData of the headerless file.
mkdir(headerFolder);
listing = dir(plainFolder);
listing = listing(~[listing.isdir]);
for k = 1:numel(listing)
    src = fullfile(plainFolder, listing(k).name);
    dst = fullfile(headerFolder, listing(k).name);
    [~, base, ext] = fileparts(listing(k).name);
    if strcmpi(ext, '.mat') && isfile(fullfile(plainFolder, [base '.dat']))
        continue   % legacy sidecar of a .dat: the headered copy has none
    end
    if ~strcmpi(ext, '.dat')
        copyfile(src, dst);
        continue
    end
    info = loadMetaData(src);
    values = loadData(src);
    hdr = struct('dataClass', info.dataClass, 'frameRateHz', info.frameRateHz, ...
        'exposureMsec', info.exposureMsec, 'channelName', base, ...
        'dimNames', {info.dimNames}, 'dimSizes', info.dimSizes, ...
        'writeComplete', true);
    fid = fopen(dst, 'w', 'ieee-le');
    fwrite(fid, encodeDatHeader(hdr), 'uint8');
    fwrite(fid, values, info.dataClass);
    fclose(fid);
end
end

function iWriteHeaderless(filePath, data, frameRateHz, exposureMsec)
% Headerless input with a legacy sidecar (.dat header Phase 5a; before it,
% headerless files described only by AcqInfos.mat).
writeTestDat(filePath, single(data), frameRateHz, exposureMsec, 'Format', 'legacySidecar');
end

function iWriteRetinotopyEvents(folder)
% Four directions, two repetitions each; every sweep is on for 10 frames
% and off for 5, so each direction has enough on-frames for the FFT bin
% used by genRetinotopyMaps (nSweeps + 1).
frameRate = 10;
directions = repmat(1:4, 1, 2);
onFrames = 1 + 15 * (0:numel(directions) - 1);
offFrames = onFrames + 10;
frames = reshape([onFrames; offFrames], [], 1);
timestamps = single(frames ./ frameRate);
state = logical(repmat([1; 0], numel(directions), 1));
eventID = uint16(repelem(directions, 2).');
repetitionID = uint16(repelem([ones(1, 4), 2 .* ones(1, 4)], 2).');
eventNameList = {'0', '180', '90', '270'};
selectedEvents = true(size(eventID));
baselinePeriod = single(0.2);
save(fullfile(folder, 'events.mat'), 'eventID', 'state', ...
    'timestamps', 'eventNameList', 'repetitionID', ...
    'selectedEvents', 'baselinePeriod');
end

function iCreateHemoComputeDataset(folder)
% Red at twice the rate of green and yellow, as in TestHemoCompute's
% repeated-illumination case, so the resampling path is exercised.
Ny = 8;
Nx = 8;
baseLength = 12;
baseFreq = 10;
AcqInfoStream = struct('Height', Ny, 'Width', Nx, 'Length', baseLength, ...
    'FrameRateHz', baseFreq, 'Camera_Model', 'D1024', 'MultiCam', false, ...
    'Datatype', 'single', 'ExposureMsec', 10);
AcqInfoStream.ImportedChannels = struct('DatFile', {}, 'Length', {}, ...
    'FrameRateHz', {}, 'ExposureMsec', {}, 'CamIdx', {});
channels = {'red.dat', 2 * baseLength, 2 * baseFreq, 0.010; ...
    'green.dat', baseLength, baseFreq, 0.015; ...
    'yellow.dat', baseLength, baseFreq, 0.020};
baseMap = single(reshape(linspace(0, 1, Ny * Nx), Ny, Nx));
for k = 1:size(channels, 1)
    nFrames = channels{k, 2};
    t = reshape(single(1:nFrames), 1, 1, []);
    data = single(1 + channels{k, 4} * sin(2*pi*t/nFrames) + 0.001 * baseMap);
    iWriteHeaderless(fullfile(folder, channels{k, 1}), data, channels{k, 3}, 10);
    AcqInfoStream.ImportedChannels(end+1) = struct('DatFile', channels{k, 1}, ...
        'Length', nFrames, 'FrameRateHz', channels{k, 3}, 'ExposureMsec', 10, 'CamIdx', 1);
end
save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
end

function store = iCreateConfiguredRig()
suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
cameras = struct('index', 1, 'displayName', 'Test Camera', 'manufacturer', '', ...
    'model', 'D1024', 'serialNumber', '', 'spectrumID', 'PF1024');
illuminations = struct('name', {'red', 'green', 'yellow'}, ...
    'displayName', {'Red', 'Green', 'Yellow'}, ...
    'manufacturer', {'', '', ''}, 'model', {'', '', ''}, ...
    'spectrumID', {'LED_632nm', 'LED_521nm', 'LED_593nm'});
store = UMITRigStore.create(struct('rigID', ['HdrEquiv_' suffix(1:10)], ...
    'cameras', cameras, 'illuminations', illuminations));
end

function iBindDatasetToRig(folder, store)
loaded = load(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
AcqInfoStream = loaded.AcqInfoStream;
rigInfo = store.getRigInfo();
AcqInfoStream.rigUUID = rigInfo.uuid;
AcqInfoStream.rigID = rigInfo.rigID;
save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream', '-mat');
end

function iRemoveFolder(folder)
if isfolder(folder)
    rmdir(folder, 's');
end
end

function out = iHemoComputeOutputs(folder)
[hbo, hbr] = HemoCompute(folder, folder, 'none', {'red', 'green', 'yellow'}, true, ...
    'RAMSafeMode', true);
out = {loadData(fullfile(folder, hbo)), loadData(fullfile(folder, hbr))};
end

function outFile = iAnaSpeckleStandardSaved(folder)
% Standard mode without outputs saves Flow.dat through saveData.
datFile = fullfile(folder, 'speckle.dat');
md = loadMetaData(datFile);
Ana_Speckle(loadData(datFile), folder, false, ...
    'FrameRateHz', double(md.frameRateHz), 'ExposureMsec', double(md.exposureMsec));
outFile = fullfile(folder, 'Flow.dat');
end

function outFile = iHemoComputeHbOPath(folder)
hbo = HemoCompute(folder, folder, 'none', {'red', 'green', 'yellow'}, true, ...
    'RAMSafeMode', true);
outFile = fullfile(folder, hbo);
end

function tf = iAnyHeadered(folder)
files = dir(fullfile(folder, '*.dat'));
tf = any(arrayfun(@(d) isDatWithHeader(fullfile(folder, d.name)), files));
end

function tf = iAllHeadered(folder)
files = dir(fullfile(folder, '*.dat'));
tf = ~isempty(files) && all(arrayfun(@(d) isDatWithHeader(fullfile(folder, d.name)), files));
end

% =========================================================================
% Output comparison
% =========================================================================
function value = iOutputValue(out, folder)
% File outputs are compared by content; in-memory outputs as returned.
value = out;
if (ischar(out) || (isstring(out) && isscalar(out))) && endsWith(out, {'.dat', '.umt'})
    filePath = char(out);
    if ~isfile(filePath)
        filePath = fullfile(folder, filePath);
    end
    value = loadData(filePath);
end
end

function filePath = iOutputPath(out, folder)
filePath = char(out);
if ~isfile(filePath)
    filePath = fullfile(folder, filePath);
end
end

function x = iNormalize(x, folder)
% Replace folder-specific text so outputs from two folders are comparable.
if ischar(x)
    x = strrep(x, folder, '<folder>');
elseif isstring(x)
    x = replace(x, folder, '<folder>');
elseif iscell(x)
    for k = 1:numel(x)
        x{k} = iNormalize(x{k}, folder);
    end
elseif isstruct(x)
    for k = 1:numel(x)
        names = fieldnames(x);
        for n = 1:numel(names)
            x(k).(names{n}) = iNormalize(x(k).(names{n}), folder);
        end
    end
elseif isdatetime(x)
    x = NaT(size(x));
end
end
