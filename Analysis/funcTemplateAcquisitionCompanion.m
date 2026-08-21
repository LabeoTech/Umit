function outFiles = funcTemplateAcquisitionCompanion(RawFolder, SaveFolder)
%FUNCTEMPLATEACQUISITIONCOMPANION Fresh SaveFolder companion template.
%
%   OUTFILES = FUNCTEMPLATEACQUISITIONCOMPANION(RAWFOLDER, SAVEFOLDER)
%   illustrates a root side-artifact importer that may accompany exactly
%   one acquisition initializer when SAVEFOLDER is fresh. It reads current
%   AcqInfos.mat and creates funcTemplateCompanion.mat without modifying
%   acquisition metadata.
%
%   INFO = FUNCTEMPLATEACQUISITIONCOMPANION('pipelineInfo') returns the
%   PipelineManager discovery contract.
%
% Files Created
%   funcTemplateCompanion.mat - Example side artifact owned by this
%                               companion function.
%
% Fresh SaveFolder Contract
%   freshSaveFolderRole='acquisition-companion' does not grant permission
%   to initialize a folder by itself. On a fresh SaveFolder, PipelineManager
%   runs the companion only after the single initializer succeeds and
%   creates valid current-schema AcqInfos.mat. A companion must never create
%   or repair AcqInfos.mat.
%
%   Copy this file to Analysis/<Category>/<YourCompanion>.m, then replace
%   the example artifact with events or another importer-owned sidecar.
%   Rename the function, pipelineInfo name, output filename, and error IDs.
%
% See also PIPELINEMANAGER, FUNCTEMPLATEACQUISITIONINITIALIZER.

defaultOutput = {'funcTemplateCompanion.mat'};

% KEEP: PipelineManager discovers modern functions through this query.
if nargin == 1 && localIsPipelineInfoQuery(RawFolder)
    outFiles = localGetPipelineInfo(defaultOutput);
    return
end

localValidateFolder(RawFolder, ...
    'Umitoolbox:funcTemplateAcquisitionCompanion:InvalidRawFolder', ...
    'RAWFOLDER');
localValidateFolder(SaveFolder, ...
    'Umitoolbox:funcTemplateAcquisitionCompanion:InvalidSaveFolder', ...
    'SAVEFOLDER');

RawFolder = char(string(RawFolder));
SaveFolder = char(string(SaveFolder));
[isLegacy, legacyMessage] = isLegacySchemaFolder(SaveFolder);
if isLegacy
    error('Umitoolbox:funcTemplateAcquisitionCompanion:MissingAcquisition', ...
        'The companion requires valid current acquisition metadata. %s', ...
        legacyMessage);
end

loaded = load(fullfile(SaveFolder, 'AcqInfos.mat'), 'AcqInfoStream');
companionArtifact = struct();
companionArtifact.SourceRawFolder = RawFolder;
companionArtifact.ImportedChannelCount = ...
    numel(loaded.AcqInfoStream.ImportedChannels);

% EDIT: Read the real auxiliary raw input and create only the artifact
% owned by this companion. Never write AcqInfos.mat here.
notesPath = fullfile(RawFolder, 'funcTemplateCompanion.txt');
if isfile(notesPath)
    companionArtifact.RawNotes = fileread(notesPath);
else
    companionArtifact.RawNotes = '';
end

save(fullfile(SaveFolder, defaultOutput{1}), ...
    'companionArtifact', '-mat');
outFiles = defaultOutput;

end

function info = localGetPipelineInfo(defaultOutput)
%LOCALGETPIPELINEINFO Declare the fresh acquisition companion contract.

info = PipelineManager.createPipelineInfo( ...
    'funcTemplateAcquisitionCompanion', ...
    'Create an acquisition side artifact after initialization.');
info.freshSaveFolderRole = 'acquisition-companion';

info = PipelineManager.addInput(info, 'RawFolder', 'RawFolder', ...
    'Folder containing auxiliary raw acquisition inputs.', ...
    'kind', 'input', 'position', 1, 'callType', 'positional', ...
    'isData', false);
info = PipelineManager.addInput(info, 'SaveFolder', 'SaveFolder', ...
    'Dataset folder initialized by the acquisition importer.', ...
    'kind', 'input', 'position', 2, 'callType', 'positional', ...
    'isData', false);

info = PipelineManager.addOutput(info, 'outFiles', 'UnknownDataType', ...
    'file', 'Side artifact created in SaveFolder.', ...
    defaultOutput, 1, 'isData', false, 'saveFileName', '');

info.notes = { ...
    'A companion cannot initialize a fresh SaveFolder by itself.'; ...
    'It runs only after the initializer creates valid AcqInfos.mat.'; ...
    'The companion must not create or repair acquisition metadata.'};

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
