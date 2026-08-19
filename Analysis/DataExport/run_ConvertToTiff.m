function outFile = run_ConvertToTiff(data, SaveFolder)
%RUN_CONVERTTOTIFF Export image-backed data to TIFF file(s).
%
%   outFile = run_ConvertToTiff(data, SaveFolder)
%   info    = run_ConvertToTiff('pipelineInfo')
%
%   Supported inputs:
%       1) Numeric Y x X x T array
%       2) Raw .dat filename storing continuous Y x X x T data
%       3) Image UMT struct
%       4) .umt filename containing one image UMT struct
%
%   Behavior:
%       - Continuous image data produce one TIFF file.
%       - Event-split image UMT data with dimensions {'Y','X','T','E'}
%         produce one TIFF file per event instance.
%       - For event-split data, a companion text file named
%         <baseName>_info.txt is created using CSV formatting. It lists the
%         generated TIFF file name, condition name, and repetition index.
%
%   Output:
%       outFile - File manifest cell array containing the generated file
%                 name(s) saved in SaveFolder. Names are derived from the
%                 input, so they are not fixed:
%                     array input        -> img_out.tif
%                     file input         -> img_<inputStem>.tif
%                     event-split input  -> one <baseName>_C<id>_R<rep>.tif
%                                           per event
%                 A <baseName>_info.txt is always written alongside.


default_Output = {'img_out.tif','img_out_info.txt'}; %#ok<NASGU>

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outFile = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'data');
addRequired(p, 'SaveFolder', @isfolder);
parse(p, data, SaveFolder);

SaveFolder = p.Results.SaveFolder;
outFile = {};

[payload, dimNames, eventInfo, baseName] = iResolveExportInput(data, SaveFolder);

assert(ismember(numel(dimNames), [3 4]), ...
    'Umitoolbox:run_ConvertToTiff:InvalidDims', ...
    'Input data must resolve to YXT or YXTE image data.');

if isequal(dimNames, {'Y','X','T'})
    tifName = [baseName '.tif'];
    iWriteTiffStack(fullfile(SaveFolder, tifName), payload);
    outFile = {tifName};
    return
end

assert(isequal(dimNames, {'Y','X','T','E'}), ...
    'Umitoolbox:run_ConvertToTiff:InvalidDims', ...
    'Event-split export requires dimensions {''Y'',''X'',''T'',''E''}.');
assert(~isempty(fieldnames(eventInfo)), ...
    'Umitoolbox:run_ConvertToTiff:MissingEventInfo', ...
    'Event-split image export requires top-level eventInfo.');
assert(isfield(eventInfo, 'eventName') && isfield(eventInfo, 'repetitionIndex') && isfield(eventInfo, 'eventID'), ...
    'Umitoolbox:run_ConvertToTiff:InvalidEventInfo', ...
    'eventInfo must contain eventID, eventName, and repetitionIndex.');

nE = size(payload, 4);
assert(numel(eventInfo.eventID) == nE && numel(eventInfo.repetitionIndex) == nE && numel(eventInfo.eventName) == nE, ...
    'Umitoolbox:run_ConvertToTiff:InvalidEventInfoLength', ...
    'eventInfo length must match the E dimension length.');

infoRows = cell(nE, 3);
for iE = 1:nE
    tifName = sprintf('%s_C%d_R%d.tif', baseName, eventInfo.eventID(iE), eventInfo.repetitionIndex(iE));
    iWriteTiffStack(fullfile(SaveFolder, tifName), payload(:,:,:,iE));
    outFile{end+1} = tifName; %#ok<AGROW>
    infoRows{iE,1} = tifName;
    infoRows{iE,2} = char(string(eventInfo.eventName{iE}));
    infoRows{iE,3} = eventInfo.repetitionIndex(iE);
end

infoName = [baseName '_info.txt'];
iWriteEventInfoText(fullfile(SaveFolder, infoName), infoRows);
outFile{end+1} = infoName;

    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            'Export image-backed data to TIFF file(s).');
        info.version = '1.0.0';

        info = PipelineManager.addInput(info, 'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            'Image-backed input to export as TIFF.', ...
            'kind', 'input', 'position', 1, 'callType', 'positional', ...
            'isData', true, 'supportsFile', true, 'dataMode', 'either');

        info = PipelineManager.addInput(info, 'SaveFolder', 'SaveFolder', ...
            'Folder where TIFF file(s) will be saved.', ...
            'kind', 'input', 'position', 2, 'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addOutput(info, 'outFile', 'ImageTimeSeries', 'file', ...
            ['Generated TIFF file manifest saved in SaveFolder. The base name is ' ...
             '''img_out'' for array input and ''img_<inputStem>'' for file input; ' ...
             'event-split input writes one ''<baseName>_C<eventID>_R<repetition>.tif'' ' ...
             'per event instead of a single ''<baseName>.tif''. The declared names ' ...
             'are the array-input, non-event case. Read the returned manifest for ' ...
             'the names actually written.'], ...
            default_Output, 1, 'isData', true, 'saveFileName', '');
    end
end

function [payload, dimNames, eventInfo, baseName] = iResolveExportInput(data, SaveFolder)
%IRESOLVEEXPORTINPUT Resolve supported inputs to export payload.

eventInfo = struct();
baseName = 'img_out';

if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty','3d'}, mfilename, 'data');
    payload = single(data);
    dimNames = {'Y','X','T'};
    return
end

if ischar(data) || (isstring(data) && isscalar(data))
    dataFile = char(string(data));
    if ~isfile(dataFile)
        altPath = fullfile(SaveFolder, dataFile);
        if isfile(altPath)
            dataFile = altPath;
        else
            error('Umitoolbox:run_ConvertToTiff:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end

    [~, stem, ext] = fileparts(dataFile);
    if ~isempty(stem)
        baseName = ['img_' stem];
    end

    switch lower(ext)
        case '.dat'
            payload = single(loadData(dataFile));
            dimNames = {'Y','X','T'};
            return
        case '.umt'
            data = loadData(dataFile);
        otherwise
            error('Umitoolbox:run_ConvertToTiff:UnsupportedInputFile', ...
                'Unsupported input file extension "%s".', ext);
    end
end

assert(isstruct(data) && isscalar(data), ...
    'Umitoolbox:run_ConvertToTiff:UnsupportedInputType', ...
    'Unsupported input type for run_ConvertToTiff.');
validateUMTStruct(data, 'requireEventInfo', false);
assert(strcmpi(char(string(data.kind)), 'image'), ...
    'Umitoolbox:run_ConvertToTiff:InvalidUMTKind', ...
    'Input UMT must have kind = "image".');

entryNames = fieldnames(data.data);
assert(~isempty(entryNames), ...
    'Umitoolbox:run_ConvertToTiff:EmptyUMTData', ...
    'Input UMT contains no image entries.');

entry = data.data.(entryNames{1});
payload = single(entry.value);
dimNames = cellstr(string(entry.dimNames));

if isfield(data, 'eventInfo')
    eventInfo = data.eventInfo;
end
end

function iWriteTiffStack(filePath, data)
%IWRITETIFFSTACK Write a YXT stack to TIFF.

data = single(data);
assert(ndims(data) == 3, 'Umitoolbox:run_ConvertToTiff:InvalidTiffPayload', ...
    'TIFF payload must be a 3D YXT array.');

if exist('ConvertToTiff', 'file') == 2
    [folderPath, fileName, ext] = fileparts(filePath);
    ConvertToTiff(folderPath, data, [fileName ext]);
    return
end

for iFrame = 1:size(data, 3)
    frame = data(:,:,iFrame);
    if iFrame == 1
        imwrite(frame, filePath, 'tif', 'Compression', 'none');
    else
        imwrite(frame, filePath, 'tif', 'WriteMode', 'append', 'Compression', 'none');
    end
end
end

function iWriteEventInfoText(filePath, rows)
%IWRITEEVENTINFOTEXT Write event export manifest in CSV-formatted text.

fid = fopen(filePath, 'w');
assert(fid ~= -1, 'Umitoolbox:run_ConvertToTiff:FileOpenFailed', ...
    'Failed to create "%s".', filePath);
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'tiffFile,conditionName,repetitionIndex\n');
for iRow = 1:size(rows,1)
    fprintf(fid, '%s,%s,%d\n', rows{iRow,1}, rows{iRow,2}, rows{iRow,3});
end
end
