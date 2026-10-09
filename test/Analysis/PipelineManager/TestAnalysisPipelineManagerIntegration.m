classdef TestAnalysisPipelineManagerIntegration < matlab.unittest.TestCase
    %TESTANALYSISPIPELINEMANAGERINTEGRATION PM execution of Analysis nodes.
    %
    % Direct Analysis unit tests do not exercise PipelineManager's source
    % resolution, RAM/file dispatch, positional argument marshalling, or leaf
    % persistence. This suite uses small synthetic fixtures and one fresh
    % SaveFolder per RAM scenario to cover those integration contracts.

    properties
        TempRoot char
    end

    methods (TestMethodSetup)
        function createTempRoot(testCase)
            fx = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TempRoot = fx.Folder;
        end
    end

    methods (Test)
        function testPMNormalizeZScoreAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('normalizeZScore', ...
                @() iCreateImageFixture(testCase.TempRoot), 'input.dat');
        end

        function testPMSpatialGaussFiltAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('spatialGaussFilt', ...
                @() iCreateImageFixture(testCase.TempRoot), 'input.dat');
        end

        function testPMApplyDetrendAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('apply_detrend', ...
                @() iCreateImageFixture(testCase.TempRoot), 'input.dat');
        end

        function testPMNormalizeLPFAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('normalizeLPF', ...
                @() iCreateImageFixture(testCase.TempRoot, 'NumFrames', 80), ...
                'input.dat');
        end

        function testPMNormalizeBSLNAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('normalizeBSLN', ...
                @() iCreateImageFixture(testCase.TempRoot), 'input.dat');
        end

        function testPMAggregateFunctionAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('apply_aggregate_function', ...
                @() iCreateImageFixture(testCase.TempRoot), 'input.dat');
        end

        function testPMGenAmplitudeMapsAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('genAmplitudeMaps', ...
                @() iCreateEventSplitFixture(testCase.TempRoot), 'input.dat', ...
                'Extensions', {'.dat'});
        end

        function testPMGenRetinotopyMapsAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('genRetinotopyMaps', ...
                @() iCreateRetinotopyFixture(testCase.TempRoot), 'input.dat');
        end

        function testPMRunConvertToTiffAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('run_ConvertToTiff', ...
                @() iCreateImageFixture(testCase.TempRoot), 'input.dat', ...
                'Extensions', {'.tif','.txt'});
        end

        function testPMRunAnaSpeckleAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('run_Ana_Speckle', ...
                @() iCreateImageFixture(testCase.TempRoot, ...
                'InputName', 'speckle.dat', 'PositiveData', true, ...
                'ExposureSpeckleMsec', 5), 'speckle.dat', ...
                'Extensions', {'.dat'});
        end

        function testPMSplitDataByEventAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('split_data_by_event', ...
                @() iCreateImageFixture(testCase.TempRoot, ...
                'NumFrames', 80, 'WithEvents', true), 'input.dat');
        end

        function testPMGenCorrelationMatrixAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('genCorrelationMatrix', ...
                @() iCreateROIFixture(testCase.TempRoot), 'input.dat');
        end

        function testPMGetDataFromROIAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('getDataFromROI', ...
                @() iCreateROIFixture(testCase.TempRoot), 'input.dat');
        end

        function testPMCalculateResponseFeaturesAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('calculateResponseFeatures', ...
                @() iCreateResponseFeatureFixture(testCase.TempRoot), ...
                'responseInput.umt', 'Extensions', {'.umt'});
        end

        function testPMGenVSMAllRamScenarios(testCase)
            testCase.verifyFileCapableNode('genVSM', ...
                @() iCreateVSMFixture(testCase.TempRoot), 'retinotopy.umt', ...
                'Extensions', {'.umt', '.roi'});
        end

        function testPMGenImageTimeSeriesUMTRamOnlyScenarios(testCase)
            testCase.verifyRamOnlyNode('genImageTimeSeriesUMT', ...
                @() iCreateImageFixture(testCase.TempRoot), 'input.dat');
        end
    end

    methods (Access = private)
        function verifyFileCapableNode(testCase, funcName, prepareFcn, inputName, varargin)
            p = inputParser;
            addParameter(p, 'Extensions', {'.dat','.umt'});
            addParameter(p, 'Parameters', struct());
            parse(p, varargin{:});

            outputs = pmCollectScenarioOutputs(prepareFcn, funcName, ...
                'Input', inputName, ...
                'Extensions', p.Results.Extensions, ...
                'PMOptions', {'Parameters', p.Results.Parameters});

            for iScenario = 1:numel(outputs)
                testCase.verifyNotEmpty(outputs(iScenario).files, ...
                    sprintf('Scenario "%s" produced no tracked output for %s.', ...
                    outputs(iScenario).scenario, funcName));
            end

            reference = outputs(1).files;
            for iScenario = 2:numel(outputs)
                testCase.verifyEqual(outputs(iScenario).files, reference, ...
                    sprintf(['%s produced %s under "%s" but %s under "%s". ' ...
                    'RAM selection must not change a node''s output identity.'], ...
                    funcName, pmFormatFileList(outputs(iScenario).files), ...
                    outputs(iScenario).scenario, pmFormatFileList(reference), ...
                    outputs(1).scenario));
            end
        end

        function verifyRamOnlyNode(testCase, funcName, prepareFcn, inputName, varargin)
            p = inputParser;
            addParameter(p, 'ExpectedFiles', {});
            parse(p, varargin{:});

            autoFolder = prepareFcn();
            before = iListDataOutputs(autoFolder);
            pm = buildPMForScenario(autoFolder, funcName, 'auto', ...
                'Input', inputName);
            result = pm.executePipeline('PrintSummary', false);
            testCase.verifyEqual(result.status, "completed", ...
                sprintf('%s failed under auto.', funcName));
            testCase.verifyNotEmpty(setdiff(iListDataOutputs(autoFolder), before));
            iVerifyExpectedFiles(testCase, result, autoFolder, ...
                p.Results.ExpectedFiles, funcName, 'auto');

            strictFolder = prepareFcn();
            pm = buildPMForScenario(strictFolder, funcName, ...
                'ramsafe-strict', 'Input', inputName);
            testCase.verifyError( ...
                @() pm.executePipeline('PrintSummary', false), ...
                'PipelineManager:RAMSafe:UnsupportedFileInput');

            bestEffortFolder = prepareFcn();
            before = iListDataOutputs(bestEffortFolder);
            pm = buildPMForScenario(bestEffortFolder, funcName, ...
                'ramsafe-bestEffort', 'Input', inputName);
            testCase.verifyWarning( ...
                @() pm.executePipeline('PrintSummary', false), ...
                'PipelineManager:RAMSafe:UnsupportedFileInput');
            testCase.verifyEqual(pm.lastExecutionResult.status, "completed", ...
                sprintf('%s failed under ramsafe-bestEffort.', funcName));
            testCase.verifyNotEmpty(setdiff( ...
                iListDataOutputs(bestEffortFolder), before));
            iVerifyExpectedFiles(testCase, pm.lastExecutionResult, ...
                bestEffortFolder, p.Results.ExpectedFiles, funcName, ...
                'ramsafe-bestEffort');
        end
    end
end

function iVerifyExpectedFiles(testCase, result, folder, expectedFiles, funcName, scenario)
for iFile = 1:numel(expectedFiles)
    fileName = string(expectedFiles{iFile});
    testCase.verifyTrue(isfile(fullfile(folder, fileName)), ...
        sprintf('%s did not create %s under %s.', funcName, fileName, scenario));
    testCase.verifyTrue(any(strcmpi(result.createdFiles.FileName, fileName)), ...
        sprintf('%s did not register %s under %s.', funcName, fileName, scenario));
end
end

function folder = iCreateImageFixture(rootFolder, varargin)
p = inputParser;
addParameter(p, 'InputName', 'input.dat');
addParameter(p, 'NumFrames', 40);
addParameter(p, 'WithEvents', false);
addParameter(p, 'PositiveData', false);
addParameter(p, 'ExposureSpeckleMsec', []);
parse(p, varargin{:});

folder = tempname(rootFolder);
mkdir(folder);

Ny = 8;
Nx = 7;
Nt = p.Results.NumFrames;
frameRate = 10;
rng(17);
data = randn(Ny, Nx, Nt, 'single');
if p.Results.PositiveData
    data = 1 + 0.05 .* data;
end
exposureMsec = 5;
if ~isempty(p.Results.ExposureSpeckleMsec)
    exposureMsec = p.Results.ExposureSpeckleMsec;
end
iWriteDat(fullfile(folder, p.Results.InputName), data, frameRate, exposureMsec);

AcqInfoStream = struct( ...
    'Height', Ny, ...
    'Width', Nx, ...
    'Length', Nt, ...
    'FrameRateHz', frameRate, ...
    'Freq', frameRate, ...
    'Datatype', 'single', ...
    'ExposureMsec', 5, ...
    'Camera_Model', 'SyntheticPMTest', ...
    'MultiCam', false, ...
    'AISampleRate', 1000, ...
    'AINChannels', 1, ...
    'AICh1', 'CameraTrig');
if ~isempty(p.Results.ExposureSpeckleMsec)
    AcqInfoStream.ExposureSpeckleMsec = p.Results.ExposureSpeckleMsec;
end
save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
ensurePMReadyAcqInfos(folder, p.Results.InputName);

if p.Results.WithEvents
    iWriteEvents(folder, Nt, frameRate);
end
end

function folder = iCreateEventSplitFixture(rootFolder)
%ICREATEEVENTSPLITFIXTURE Folder whose input.dat is event-split (Y-X-T-E).
folder = iCreateImageFixture(rootFolder, 'NumFrames', 80, 'WithEvents', true);
splitFile = char(string(split_data_by_event(fullfile(folder, 'input.dat'), folder)));
movefile(splitFile, fullfile(folder, 'input.dat'), 'f');
end

function folder = iCreateRetinotopyFixture(rootFolder)
folder = iCreateImageFixture(rootFolder, 'NumFrames', 64);
eventID = uint16([1;1;2;2;3;3;4;4;1;1;2;2;3;3;4;4]);
state = logical(repmat([1;0], 8, 1));
% Eight 0.4 s sweeps (4 frames at 10 Hz), 0.5 s apart, leaving room before
% the first onset for the average-movie baseline.
onsetSec = 1.8 + 0.5 .* (0:7).';
timestamps = single(reshape([onsetSec, onsetSec + 0.4].', [], 1));
eventNameList = {'0','180','90','270'};
repetitionID = uint16([ones(8,1); 2.*ones(8,1)]);
selectedEvents = true(size(eventID));
baselinePeriod = single(0.1);
save(fullfile(folder, 'events.mat'), 'eventID', 'state', ...
    'timestamps', 'eventNameList', 'repetitionID', ...
    'selectedEvents', 'baselinePeriod');
end

function folder = iCreateROIFixture(rootFolder)
folder = iCreateImageFixture(rootFolder);
imageSizeYX = [8 7];
maskA = false(imageSizeYX); maskA(2,2) = true;
maskB = false(imageSizeYX); maskB(5,5) = true;
rois = [iMakeROI('ROI_A', maskA), iMakeROI('ROI_B', maskB)];
ROIFile = createROIFile(imageSizeYX, 'ROIs', rois);
saveROIFile(fullfile(folder, 'myROI.roi'), ROIFile);
end

function folder = iCreateResponseFeatureFixture(rootFolder)
folder = iCreateImageFixture(rootFolder);
vals = zeros(2, 12, 2, 'single');
vals(1,:,1) = [0 0 0 0 0 1 2 3 2 1 0 0];
vals(1,:,2) = [0 0 0 0 0 2 3 4 3 2 0 0];
vals(2,:,:) = 0.5 .* vals(1,:,:);
labels = struct('ROI', {{'ROI1','ROI2'}});
umt = genUMTStruct(vals, 'kind', 'roi', 'entryName', 'main', ...
    'dimNames', {'ROI','T','E'}, 'labels', labels);
umt.data.main.meta = struct('FrameRateHz', 10);
umt = appendUMTEventInfo(umt, ...
    'eventID', [1;2], 'repetitionIndex', [1;1], ...
    'eventName', {'A';'B'}, 'eventAxisMode', 'instances', ...
    'overwrite', true);
umt.eventInfo.baselinePeriod = 0.5;
saveData(fullfile(folder, 'responseInput.umt'), umt);
end

function folder = iCreateVSMFixture(rootFolder)
folder = iCreateImageFixture(rootFolder);
[y, x] = ndgrid(single(1:8), single(1:7));
az = cat(3, ones(8,7,'single'), x ./ 7);
el = cat(3, ones(8,7,'single'), y ./ 8);
umt = genUMTStruct(az, 'kind', 'image', 'entryName', ...
    'AzimuthMap', 'dimNames', {'Y','X','F'});
umt = genUMTStruct(umt, 'value', el, 'entryName', ...
    'ElevationMap', 'dimNames', {'Y','X','F'});
saveData(fullfile(folder, 'retinotopy.umt'), umt);
end

function roi = iMakeROI(name, mask)
[rows, cols] = find(mask);
x0 = min(cols) - 0.5; x1 = max(cols) + 0.5;
y0 = min(rows) - 0.5; y1 = max(rows) + 0.5;
roi = struct();
roi.name = name;
vertices = [x0 y0; x1 y0; x1 y1; x0 y1];
pgon = polyshape(vertices);
roi.type = 'polygon';
roi.DOC = datetime('now');
roi.modifiedOn = datetime('now');
roi.color = [1 0 0];
roi.notes = '';
roi.geometry = struct( ...
    'polyshape', pgon, ...
    'verticesXY_px', vertices, ...
    'ROIType', 'polygon', ...
    'ROIParameters', struct('ROIType', 'polygon', ...
        'Position', vertices, 'Vertices', vertices));
roi.mask = logical(mask);
roi.stats = struct( ...
    'computedOn', datetime('now'), ...
    'NPixels', nnz(mask), ...
    'areaPx2', nnz(mask), ...
    'areaMM2', [], ...
    'centroidXY_px', [], ...
    'centroidXY_mm', [], ...
    'distanceFromOrigin_px', [], ...
    'distanceFromOrigin_mm', [], ...
    'spatialMean', [], ...
    'spatialStd', [], ...
    'spatialMedian', [], ...
    'spatialMin', [], ...
    'spatialMax', []);
end

function iWriteEvents(folder, nFrames, frameRate)
onsets = [10 30 50 70];
onsets = onsets(onsets + 4 <= nFrames);
offsets = onsets + 4;
timestamps = reshape([onsets; offsets], [], 1);
timestamps = single((timestamps - 1) ./ frameRate);
state = logical(repmat([1;0], numel(onsets), 1));
eventID = uint16(repelem(1:max(1, numel(onsets)/2), 4).');
eventID = eventID(1:numel(timestamps));
eventNameList = arrayfun(@(x) sprintf('Event%d', x), ...
    unique(eventID), 'UniformOutput', false);
selectedEvents = true(size(eventID));
baselinePeriod = single(0.2);
save(fullfile(folder, 'events.mat'), 'timestamps', 'state', ...
    'eventID', 'eventNameList', 'baselinePeriod', 'selectedEvents');
end

function iWriteDat(filePath, data, frameRateHz, exposureMsec)
% Headered input (.dat header Phase 5a).
writeTestDat(filePath, single(data), frameRateHz, exposureMsec);
end

function files = iListDataOutputs(folder)
files = {};
for ext = {'.dat','.umt'}
    found = dir(fullfile(folder, ['*' ext{1}]));
    files = [files, {found.name}]; %#ok<AGROW>
end
files = unique(files);
end
