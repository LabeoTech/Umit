function outData = funcTemplateUMT(data, SaveFolder, varargin)
%FUNCTEMPLATEUMT Template for creating a UMT pipeline output.
%
%   OUTDATA = FUNCTEMPLATEUMT(DATA, SAVEFOLDER) calculates the temporal
%   mean of a numeric image time series and packages the resulting Y-by-X
%   image in a validated UMT structure. UMT is umIT's schema-validated
%   container for processed and derived data.
%
%   OUTDATA = FUNCTEMPLATEUMT(DATA, SAVEFOLDER, Name, Value, ...) selects
%   a custom UMT entry name.
%
%   INFO = FUNCTEMPLATEUMT('pipelineInfo') returns the metadata used by
%   PipelineManager to discover and run the function.
%
% Inputs
%   data       - Nonempty, real numeric Y-by-X-by-T array held in RAM.
%   SaveFolder - Existing save-folder path supplied by PipelineManager.
%                genUMTStruct uses it to resolve applicable entry metadata,
%                such as FrameRateHz, from authoritative local files.
%
% Name-Value Option
%   EntryName - Valid MATLAB field name used for the output UMT entry.
%               Default: 'MeanImage'.
%
% Output
%   outData - Valid image UMT containing one Y-by-X entry. PipelineManager
%             passes the structure downstream and saves it as
%             "funcTemplateMean.umt" when persistence is required.
%
% Metadata and Dimension Behavior
%   The calculation removes the T dimension, so the output declares
%   dimNames={'Y','X'} instead of passing input metadata through. Supply
%   SaveFolder to genUMTStruct so supported entry metadata can be resolved
%   locally. When a copied function creates event-split data or changes
%   other dimensions, update dimNames, labels, eventInfo, and entry meta to
%   match the result before calling validateUMTStruct.
%
% Files Created or Modified
%   Direct calls create no files. This function returns a UMT structure;
%   PipelineManager manages the declared "funcTemplateMean.umt" file.
%
% Example
%   data = single(reshape(1:24, [2 3 4]));
%   saveFolder = tempdir;
%   outData = funcTemplateUMT(data, saveFolder, 'EntryName', 'Average');
%
% PipelineManager Compatibility
%   Copy this file to Analysis/<Category>/<YourFunction>.m, then rename the
%   function, pipelineInfo name, entry, and output file. Preserve the
%   one-argument 'pipelineInfo' query. This example is RAM-only and returns
%   ProcessedData through outputMode='data'.
%
% See also GENUMTSTRUCT, VALIDATEUMTSTRUCT, SAVEDATA, PIPELINEMANAGER.

% KEEP: PipelineManager discovers modern functions through this query.
if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localGetPipelineInfo();
    return
end

% EDIT: Add options here and mirror them in localGetPipelineInfo.
p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'data');
addRequired(p, 'SaveFolder');
addParameter(p, 'EntryName', 'MeanImage', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p, data, SaveFolder, varargin{:});

data = p.Results.data;
SaveFolder = p.Results.SaveFolder;
entryName = char(string(p.Results.EntryName));

localValidateData(data);
localValidateSaveFolder(SaveFolder);
if ~isvarname(entryName)
    error('umIToolbox:funcTemplateUMT:InvalidEntryName', ...
        'ENTRYNAME must be a valid MATLAB field name.');
end

% EDIT: Replace this minimal calculation with the scientific operation.
meanImage = mean(data, 3, 'omitnan');

% KEEP OR UPDATE TOGETHER: the UMT kind and dimNames must describe the
% returned value. SaveFolder enables local metadata resolution.
outData = genUMTStruct( ...
    meanImage, ...
    'kind', 'image', ...
    'entryName', entryName, ...
    'dimNames', {'Y','X'}, ...
    'SaveFolder', SaveFolder);

validateUMTStruct(outData, 'requireEventInfo', false);

end

function info = localGetPipelineInfo()
%LOCALGETPIPELINEINFO Declare the PipelineManager integration contract.

info = PipelineManager.createPipelineInfo( ...
    'funcTemplateUMT', ...
    'Calculate a temporal-mean image and return it in UMT format.');

info = PipelineManager.addInput( ...
    info, ...
    'data', ...
    'ImageTimeSeries', ...
    'Numeric YXT image time series held in RAM.', ...
    'kind', 'input', ...
    'position', 1, ...
    'callType', 'positional', ...
    'isData', true, ...
    'supportsFile', false, ...
    'dataMode', 'ram');

info = PipelineManager.addInput( ...
    info, ...
    'SaveFolder', ...
    'SaveFolder', ...
    'Folder used for local metadata resolution.', ...
    'kind', 'input', ...
    'position', 2, ...
    'callType', 'positional', ...
    'isData', false);

info = PipelineManager.addInput( ...
    info, ...
    'EntryName', ...
    'parameter', ...
    'Valid field name for the generated UMT entry.', ...
    'kind', 'parameter', ...
    'position', 3, ...
    'callType', 'namevalue', ...
    'default', 'MeanImage', ...
    'dataType', 'char');

info = PipelineManager.addOutput( ...
    info, ...
    'outData', ...
    'ProcessedData', ...
    'data', ...
    'Temporal-mean image in a validated UMT structure.', ...
    'funcTemplateMean.umt', ...
    1, ...
    'isData', true, ...
    'saveFileName', 'funcTemplateMean.umt');

info.notes = { ...
    'Use genUMTStruct for creation or appending and validate the result.'; ...
    'PipelineManager, not this function, persists the returned UMT data.'};

end

function localValidateData(data)
%LOCALVALIDATEDATA Check the documented YXT input contract.

if ~(isnumeric(data) && isreal(data) && ~isempty(data) && ndims(data) == 3)
    error('umIToolbox:funcTemplateUMT:InvalidData', ...
        'DATA must be a nonempty, real numeric Y-by-X-by-T array.');
end
end

function localValidateSaveFolder(SaveFolder)
%LOCALVALIDATESAVEFOLDER Require an existing dataset save folder.

isTextScalar = ischar(SaveFolder) || ...
    (isstring(SaveFolder) && isscalar(SaveFolder));
if ~isTextScalar || ~isfolder(SaveFolder)
    error('umIToolbox:funcTemplateUMT:InvalidSaveFolder', ...
        'SAVEFOLDER must be the path to an existing folder.');
end
end
