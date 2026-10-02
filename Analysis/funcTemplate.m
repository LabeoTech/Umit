function outData = funcTemplate(data, SaveFolder, varargin)
%FUNCTEMPLATE Template for a PipelineManager numeric analysis function.
%
%   OUTDATA = FUNCTEMPLATE(DATA, SAVEFOLDER) multiplies a numeric image
%   time series by a configurable scale factor. This deliberately simple
%   operation keeps the example focused on PipelineManager integration.
%
%   OUTDATA = FUNCTEMPLATE(DATA, SAVEFOLDER, Name, Value, ...) applies the
%   operation with user-selected options.
%
%   INFO = FUNCTEMPLATE('pipelineInfo') returns the metadata used by
%   PipelineManager to discover, validate, connect, and run the function.
%   A semantic type, such as ImageTimeSeries, describes what a pipeline
%   port carries; it is not a MATLAB class.
%
% Inputs
%   data       - Nonempty, real numeric Y-by-X-by-T array held in RAM.
%                Y and X are image dimensions and T is time. This template
%                does not accept a .dat filename.
%   SaveFolder - Existing save-folder path supplied by PipelineManager.
%                Use this path to resolve files and metadata owned by the
%                current dataset.
%
% Name-Value Options
%   Scale        - Finite, nonnegative numeric scalar multiplied into DATA.
%                  Default: 1. The declared range [0 Inf] is displayed as
%                  an editable PipelineManager option.
%   ClipNegative - Logical scalar. When true, negative output samples are
%                  replaced with zero. Default: false.
%
% Output
%   outData - Numeric Y-by-X-by-T result with the same dimensions and
%             MATLAB class as DATA. PipelineManager can pass this value to
%             a connected downstream data input.
%
% Metadata Behavior
%   Per-data metadata (frame rate, axes, exposure) describes the data, not
%   the pipeline. A function that needs it declares 'sourceInfo' inputs
%   (FrameRateHz, DimNames, ExposureMsec; see PipelineManager.addInput),
%   which PipelineManager injects as Name-Value parameters from the data
%   flowing into the step. Inside the function, resolve the value with
%   resolveDataInfoValue (explicit value, then the .dat header, then a
%   clear error). Never read the frame rate from AcqInfos.mat.
%
%   When a copied function accepts a .dat filename, describe the file and
%   check its layout before reading data or writing outputs, for example:
%
%       dataFile = fullfile(SaveFolder, dataFileName);
%       dataInfo = loadMetaData(dataFile);
%       assertDatLayout(dataInfo, {{'Y','X','T'}}, mfilename);
%       nFrames  = datAxisSize(dataInfo, 'T');   % 0 when the file has no T axis
%
%   .dat files may be Y-X, Y-X-T, Y-X-T-E, Y-X-E, or Y-X-F; a function that
%   assumes Y-X-T must refuse the others (Umitoolbox:<fn>:unsupportedLayout).
%   Outputs saved by PipelineManager take their axes from the source Info
%   or the declared output type; a direct saveData call must pass
%   'DimNames'.
%
%   For .dat files, read only the .dat Info schema fields: dimNames,
%   dimSizes, dataClass, frameRateHz, exposureMsec, dataOffset, format,
%   filePath, channelName, writeComplete. Use datAxisSize for single axis
%   sizes. The older names (Height, Width, Length, datLength, datSize,
%   dim_names, Datatype, Freq, FrameRateHz, ...) are deprecated and will be
%   removed. To stream parts of a large .dat file, use
%   spatialSlabIO('open'/'read'/'close') instead of reading the file.
%
%   Use the file associated with the data being processed; do not select an
%   arbitrary file from SaveFolder. Do not edit AcqInfos.mat or another
%   authoritative metadata file from an ordinary processing function.
%
% Files Created or Modified
%   Direct calls create no files. PipelineManager may save OUTDATA as
%   "funcTemplate.dat", the declared managed output name. Replace that
%   name when renaming a copied function.
%
% Example
%   data = single(reshape(-12:11, [2 3 4]));
%   saveFolder = tempdir;
%   outData = funcTemplate(data, saveFolder, ...
%       'Scale', 2, 'ClipNegative', true);
%
% PipelineManager Compatibility
%   Copy this file to Analysis/<Category>/<YourFunction>.m. Rename the
%   function, update the name passed to createPipelineInfo, and choose a
%   unique default output filename. Preserve the one-argument
%   'pipelineInfo' query and keep argument positions aligned with the
%   MATLAB function declaration. Declare supportsFile=true only after
%   implementing and testing filename input. This example has no external
%   dependencies.
%
% See also PIPELINEMANAGER, PIPELINEMANAGER.CREATEPIPELINEINFO,
% PIPELINEMANAGER.ADDINPUT, PIPELINEMANAGER.ADDOUTPUT, LOADMETADATA.

% KEEP: PipelineManager discovers modern functions through this query.
if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localGetPipelineInfo();
    return
end

% EDIT: Add user-facing options here and mirror them in localGetPipelineInfo.
p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'data');
addRequired(p, 'SaveFolder');
addParameter(p, 'Scale', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'ClipNegative', false, ...
    @(x) islogical(x) && isscalar(x));
parse(p, data, SaveFolder, varargin{:});

data = p.Results.data;
SaveFolder = p.Results.SaveFolder;
scale = p.Results.Scale;
clipNegative = p.Results.ClipNegative;

localValidateData(data);
localValidateSaveFolder(SaveFolder);

% EDIT: Put the scientific calculation in this section. Cast option values
% like DATA when the result should preserve its MATLAB numeric class.
outData = data .* cast(scale, 'like', data);
if clipNegative
    outData(outData < 0) = 0;
end

end

function info = localGetPipelineInfo()
%LOCALGETPIPELINEINFO Declare the PipelineManager integration contract.

% CONVENTION: new Analysis/ functions should use the positional order
% (data, SaveFolder, ...) shown by this template's own signature. A few
% existing functions predate this convention and use a different order;
% follow this template for new functions instead of copying them.
%
% CONVENTION: a parameter that accepts either a keyword or a numeric range
% (e.g. 'auto' or an explicit [min max]) should declare 'allowed' as
% {'keyword', [min max]}, matching calculateResponseFeatures'
% TimeWindow_sec and genAmplitudeMaps' equivalent parameter. Do not invent
% a new encoding for this shape.

% EDIT: Rename the function, description, and output file in every copy.
info = PipelineManager.createPipelineInfo( ...
    'funcTemplate', ...
    'Multiply an in-memory YXT image time series by a scale factor.');

% KEEP THE STRUCTURE: each positional function input needs one input entry.
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

% SaveFolder is a dedicated non-data input. PipelineManager supplies the
% current dataset folder automatically during execution.
info = PipelineManager.addInput( ...
    info, ...
    'SaveFolder', ...
    'SaveFolder', ...
    'Folder used to resolve dataset files and metadata.', ...
    'kind', 'input', ...
    'position', 2, ...
    'callType', 'positional', ...
    'isData', false);

% EDIT TOGETHER: defaults, accepted values, and validators must agree.
info = PipelineManager.addInput( ...
    info, ...
    'Scale', ...
    'parameter', ...
    'Nonnegative scale factor applied to every sample.', ...
    'kind', 'parameter', ...
    'position', 3, ...
    'callType', 'namevalue', ...
    'default', 1, ...
    'allowed', [0 Inf], ...
    'dataType', 'double');

info = PipelineManager.addInput( ...
    info, ...
    'ClipNegative', ...
    'parameter', ...
    'Replace negative output samples with zero.', ...
    'kind', 'parameter', ...
    'position', 4, ...
    'callType', 'namevalue', ...
    'default', false, ...
    'allowed', [true false], ...
    'dataType', 'logical');

% outputMode='data' means the function returns the value. PipelineManager
% passes it downstream and manages persistence when required.
info = PipelineManager.addOutput( ...
    info, ...
    'outData', ...
    'ImageTimeSeries', ...
    'data', ...
    'Scaled YXT image time series.', ...
    'funcTemplate.dat', ...
    1, ...
    'isData', true, ...
    'saveFileName', 'funcTemplate.dat');

info.notes = { ...
    'RAM-only example: supportsFile is intentionally false.'; ...
    'Resolve required metadata locally from a concrete file in SaveFolder.'};

end

function localValidateData(data)
%LOCALVALIDATEDATA Check the documented YXT input contract.

if ~(isnumeric(data) && isreal(data) && ~isempty(data) && ndims(data) == 3)
    error('Umitoolbox:funcTemplate:InvalidData', ...
        'DATA must be a nonempty, real numeric Y-by-X-by-T array.');
end
end

function localValidateSaveFolder(SaveFolder)
%LOCALVALIDATESAVEFOLDER Require an existing dataset save folder.

isTextScalar = ischar(SaveFolder) || ...
    (isstring(SaveFolder) && isscalar(SaveFolder));
if ~isTextScalar || ~isfolder(SaveFolder)
    error('Umitoolbox:funcTemplate:InvalidSaveFolder', ...
        'SAVEFOLDER must be the path to an existing folder.');
end
end
