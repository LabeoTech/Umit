function outFiles = funcTemplateAcquisitionInitializer(RawFolder, SaveFolder)
%FUNCTEMPLATEACQUISITIONINITIALIZER Fresh SaveFolder importer template.
%
%   OUTFILES = FUNCTEMPLATEACQUISITIONINITIALIZER(RAWFOLDER, SAVEFOLDER)
%   illustrates an acquisition importer that may initialize an unprocessed
%   SAVEFOLDER. RAWFOLDER and SAVEFOLDER may be the same directory. For a
%   runnable example, RAWFOLDER must contain
%   funcTemplateAcquisition.mat with variables imageData and
%   AcqInfoStream. The metadata must describe IMAGEData exactly.
%
%   INFO = FUNCTEMPLATEACQUISITIONINITIALIZER('pipelineInfo') returns the
%   PipelineManager discovery contract.
%
% Files Created
%   funcTemplateImported.dat - Imported Y-X-T image data, written with
%                              saveData(..., 'FrameRateHz', rate, 'Info',
%                              struct('exposureMsec', exposure)) so the
%                              file carries its own .dat header.
%   AcqInfos.mat             - Current acquisition metadata owned by this
%                              importer, including ImportedChannels.
%
% Fresh SaveFolder Contract
%   freshSaveFolderRole='acquisition-initializer' is reserved for an
%   importer that creates authoritative AcqInfos.mat metadata from raw
%   acquisition facts. PipelineManager does not create or repair that
%   metadata. The initializer must return successfully only after writing
%   a valid current-schema AcqInfos.mat.
%
%   Copy this file to Analysis/<Category>/<YourImporter>.m, then replace
%   the example MAT-file reader with the real acquisition reader. Rename
%   the function, pipelineInfo name, output filename, and error IDs.
%
% See also PIPELINEMANAGER, FUNCTEMPLATEACQUISITIONCOMPANION,
% APPENDIMPORTEDCHANNELINFO.

defaultOutput = {'funcTemplateImported.dat'};

% KEEP: PipelineManager discovers modern functions through this query.
if nargin == 1 && localIsPipelineInfoQuery(RawFolder)
    outFiles = localGetPipelineInfo(defaultOutput);
    return
end

localValidateFolder(RawFolder, ...
    'Umitoolbox:funcTemplateAcquisitionInitializer:InvalidRawFolder', ...
    'RAWFOLDER');
localValidateFolder(SaveFolder, ...
    'Umitoolbox:funcTemplateAcquisitionInitializer:InvalidSaveFolder', ...
    'SAVEFOLDER');
localRequireUninitializedSaveFolder(SaveFolder);

RawFolder = char(string(RawFolder));
SaveFolder = char(string(SaveFolder));
sourcePath = fullfile(RawFolder, 'funcTemplateAcquisition.mat');
if ~isfile(sourcePath)
    error('Umitoolbox:funcTemplateAcquisitionInitializer:MissingSource', ...
        ['The example source file "%s" is missing. Replace this reader ' ...
         'with the real acquisition-format reader.'], sourcePath);
end

source = load(sourcePath, 'imageData', 'AcqInfoStream');
if ~isfield(source, 'imageData') || ~isfield(source, 'AcqInfoStream')
    error('Umitoolbox:funcTemplateAcquisitionInitializer:InvalidSource', ...
        'The example source must contain imageData and AcqInfoStream.');
end

imageData = source.imageData;
AcqInfoStream = source.AcqInfoStream;
localValidateAcquisition(imageData, AcqInfoStream);

% EDIT: Import the real raw data before publishing its metadata.
dataPath = fullfile(SaveFolder, defaultOutput{1});
localWriteImageData(dataPath, imageData, AcqInfoStream);

% Importers own AcqInfos.mat. Build ImportedChannels from the raw facts;
% do not ask PipelineManager or a companion function to synthesize it.
if isfield(AcqInfoStream, 'ImportedChannels')
    AcqInfoStream = rmfield(AcqInfoStream, 'ImportedChannels');
end
channelInfo = struct( ...
    'DatFile', defaultOutput{1}, ...
    'Length', size(imageData, 3), ...
    'FrameRateHz', AcqInfoStream.FrameRateHz);
if isfield(AcqInfoStream, 'ExposureMsec')
    channelInfo.ExposureMsec = AcqInfoStream.ExposureMsec;
end
AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, channelInfo);
save(fullfile(SaveFolder, 'AcqInfos.mat'), 'AcqInfoStream', '-mat');

outFiles = defaultOutput;

end

function info = localGetPipelineInfo(defaultOutput)
%LOCALGETPIPELINEINFO Declare the fresh acquisition initializer contract.

info = PipelineManager.createPipelineInfo( ...
    'funcTemplateAcquisitionInitializer', ...
    'Import acquisition data and initialize a fresh SaveFolder.');
info.freshSaveFolderRole = 'acquisition-initializer';

info = PipelineManager.addInput(info, 'RawFolder', 'RawFolder', ...
    'Folder containing authoritative raw acquisition inputs.', ...
    'kind', 'input', 'position', 1, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addInput(info, 'SaveFolder', 'SaveFolder', ...
    'Unprocessed destination initialized by this importer.', ...
    'kind', 'input', 'position', 2, 'callType', 'positional', ...
    'isData', false);

info = PipelineManager.addOutput(info, 'outFiles', 'ImageTimeSeries', ...
    'file', 'Imported image files created in SaveFolder.', ...
    defaultOutput, 1, 'isData', true, 'saveFileName', '');

info.notes = { ...
    'Exactly one acquisition initializer is allowed for a fresh SaveFolder.'; ...
    'The function must create valid current-schema AcqInfos.mat metadata.'; ...
    'PipelineManager never creates or repairs acquisition metadata.'};

end

function tf = localIsPipelineInfoQuery(value)
%LOCALISPIPELINEINFOQUERY Recognize the discovery-only call.

tf = (ischar(value) || (isstring(value) && isscalar(value))) && ...
    strcmpi(strtrim(char(string(value))), 'pipelineInfo');
end

function localValidateFolder(value, errorID, argumentName)
%LOCALVALIDATEFOLDER Require an existing folder path.

isTextScalar = ischar(value) || (isstring(value) && isscalar(value));
if ~isTextScalar || ~isfolder(value)
    error(errorID, '%s must be the path to an existing folder.', argumentName);
end
end

function localRequireUninitializedSaveFolder(SaveFolder)
%LOCALREQUIREUNINITIALIZEDSAVEFOLDER Avoid overwriting dataset artifacts.

blockedExtensions = {'.dat', '.umt', '.roi', '.umitlink'};
blockedNames = { ...
    'AcqInfos.mat', 'DataParams.mat', 'events.mat', 'Text_events.mat', ...
    'dataHistory.mat', 'pipeLog.mat', 'LogBook.mat'};
entries = dir(SaveFolder);
for iEntry = 1:numel(entries)
    if entries(iEntry).isdir
        continue
    end

    entryName = entries(iEntry).name;
    [~, ~, extension] = fileparts(entryName);
    if any(strcmpi(entryName, blockedNames)) || ...
            any(strcmpi(extension, blockedExtensions))
        error('Umitoolbox:funcTemplateAcquisitionInitializer:InitializedSaveFolder', ...
            ['This example importer requires a SAVEFOLDER with no existing ' ...
             'processed-dataset artifacts. Raw acquisition files are allowed.']);
    end
end
end

function localValidateAcquisition(imageData, AcqInfoStream)
%LOCALVALIDATEACQUISITION Check that raw facts describe the imported data.

if ~(isnumeric(imageData) && isa(imageData, 'single') && ...
        isreal(imageData) && ~isempty(imageData) && ...
        ndims(imageData) == 3)
    error('Umitoolbox:funcTemplateAcquisitionInitializer:InvalidImageData', ...
        'imageData must be a nonempty, real single Y-X-T array.');
end
if ~(isstruct(AcqInfoStream) && isscalar(AcqInfoStream))
    error('Umitoolbox:funcTemplateAcquisitionInitializer:InvalidMetadata', ...
        'AcqInfoStream must be a scalar structure.');
end

requiredFields = {'Width', 'Height', 'Length', 'FrameRateHz'};
if ~all(isfield(AcqInfoStream, requiredFields))
    error('Umitoolbox:funcTemplateAcquisitionInitializer:InvalidMetadata', ...
        'AcqInfoStream must contain Width, Height, Length, and FrameRateHz.');
end
expectedSize = [double(AcqInfoStream.Height), ...
    double(AcqInfoStream.Width), double(AcqInfoStream.Length)];
if ~isequal(double(size(imageData)), expectedSize) || ...
        ~(isnumeric(AcqInfoStream.FrameRateHz) && ...
        isscalar(AcqInfoStream.FrameRateHz) && ...
        isfinite(AcqInfoStream.FrameRateHz) && ...
        AcqInfoStream.FrameRateHz > 0)
    error('Umitoolbox:funcTemplateAcquisitionInitializer:InvalidMetadata', ...
        'AcqInfoStream dimensions and frame rate must describe imageData.');
end
end

function localWriteImageData(dataPath, imageData, AcqInfoStream)
%LOCALWRITEIMAGEDATA Write the imported stream before publishing metadata.
%
% Imported .dat files are written with saveData, which adds the
% self-describing header. Pass the stream's own frame rate ('FrameRateHz')
% and exposure ('Info' with exposureMsec) so the file does not depend on
% AcqInfos.mat to be read. AcqInfos.mat itself is written by the caller.

exposureMsec = NaN;
if isfield(AcqInfoStream, 'ExposureMsec') && ~isempty(AcqInfoStream.ExposureMsec)
    exposureMsec = double(AcqInfoStream.ExposureMsec);
end
saveData(dataPath, imageData, ...
    'FrameRateHz', double(AcqInfoStream.FrameRateHz), ...
    'Info', struct('exposureMsec', exposureMsec));
end
