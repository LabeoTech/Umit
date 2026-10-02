function outFile = saveData(filename, data, varargin)
%SAVEDATA Save image data or derived UMIT data to disk.
%
%   saveData(filename, data, 'DimNames', dimNames)
%   saveData(..., 'Info', Info)
%   saveData(..., 'FrameRateHz', rate)
%   saveData(..., 'ChannelName', name)
%   saveData(..., 'Append', true)
%   saveData(filename, umtStruct)
%
%   Inputs:
%       filename        - Full path of the file to be saved. The extension
%                         is forced based on the input data type.
%       data            - Either:
%                         1) non-empty single numeric array, saved as a
%                            headered .dat file with the axes in DimNames
%                         2) valid derived-data structure, saved as .umt
%                            (its entries carry their own dimNames)
%
%   Name-Value options (numeric data only):
%       DimNames        - Required. Axis names of DATA in memory order: 'Y',
%                         'X', then any of 'T', 'E', 'F' in that order, e.g.
%                         {'Y','X'}, {'Y','X','T'}, {'Y','X','T','E'},
%                         {'Y','X','E'}. saveData never guesses the layout.
%                         DATA must fit them: ndims(data) <= numel(DimNames)
%                         (MATLAB drops trailing singleton axes, so a 2-D
%                         array can be saved as Y-X-T with T = 1), and the
%                         axis sizes are size(data, 1:numel(DimNames)).
%       Info            - .dat Info struct of the source file (as returned
%                         by loadMetaData). Its frameRateHz and exposureMsec
%                         fields, when valid, describe the output.
%       FrameRateHz     - Positive finite frame rate of the output.
%                         Overrides Info.frameRateHz.
%       ChannelName     - Header channelName. Default: the output file's
%                         base name. (PipelineManager passes the save
%                         target's name for temporary files that may be
%                         promoted to it.)
%       Append          - Logical scalar. If true, append DATA to an
%                         existing headered .dat file along its last axis.
%                         Default: false
%
%   Notes:
%       - Numeric arrays are written with a version-1 header (see
%         docs/dev/dat-header-spec.md): axes DimNames, class single,
%         channelName = the output file's base name (non-ASCII characters
%         replaced by '_', truncated to the header field).
%       - Frame rate, with a T axis: FrameRateHz, else Info.frameRateHz;
%         with neither, Umitoolbox:saveData:missingFrameRate is raised and
%         nothing is written. Without a T axis the header frame rate is NaN
%         and both are ignored. AcqInfos.mat is never read or written: it
%         describes the raw acquisition, whose rate can differ from the
%         imported data (temporal binning).
%       - Exposure: Info.exposureMsec, else NaN.
%       - An existing file is overwritten, headered or not.
%       - Append to a missing file creates it. Append to a headered file
%         requires the same axes (DimNames), class, and frame rate, the
%         same size on every axis except the last, and a complete earlier
%         write; the last axis grows by DATA's size along it. Append to a
%         headerless file is an error and leaves the file unchanged.
%       - Structured derived data are saved as MAT-files with ".umt".

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'saveData';

addRequired(p, 'filename', @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(p, 'data', @(x) (isnumeric(x) && isa(x, 'single') && ~isempty(x)) || isstruct(x));
addParameter(p, 'DimNames', {}, @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(p, 'Append', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Info', struct.empty(0,1), @(x) isstruct(x) || isempty(x));
addParameter(p, 'FrameRateHz', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x > 0));
addParameter(p, 'ChannelName', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));

parse(p, filename, data, varargin{:});

filename = convertStringsToChars(filename);
[saveFolder, fileBase, ~] = fileparts(filename);

if isempty(saveFolder)
    saveFolder = pwd;
end

if isstruct(data)
    outFile = [fileBase,'.umt'];
    save2umt(fullfile(saveFolder, outFile), data);
else
    outFile = [fileBase, '.dat'];
    save2dat(fullfile(saveFolder,outFile), data, p.Results);
end

disp(['Data saved as "' outFile '"']);
end

% =========================================================================
% Local functions
% =========================================================================

function save2dat(filePath, data, opts)
%SAVE2DAT Save a numeric array to a headered .dat file with explicit axes.

if ~(isnumeric(data) && isa(data, 'single') && ~isempty(data))
    error('Umitoolbox:saveData:invalidInput', ...
        'Numeric data must be a non-empty single array.');
end

[dimNames, dimSizes] = iResolveLayout(opts.DimNames, data, filePath);

[saveFolder, fileName, ext] = fileparts(filePath);
if isempty(saveFolder)
    saveFolder = pwd;
end
if ~strcmpi(ext, '.dat')
    filePath = fullfile(saveFolder, [fileName, '.dat']);
end

if ~isfolder(saveFolder)
    error('Umitoolbox:saveData:invalidFolder', ...
        'Target folder does not exist: "%s".', saveFolder);
end

[frameRateHz, exposureMsec] = iResolveRateAndExposure(opts, filePath, any(strcmp(dimNames, 'T')));

channelName = fileName;
if strlength(string(opts.ChannelName)) > 0
    channelName = char(opts.ChannelName);
end
hdr = datHeaderFromInfo(struct('dataClass', 'single', 'dimNames', {dimNames}, ...
    'dimSizes', dimSizes, 'frameRateHz', frameRateHz, 'exposureMsec', exposureMsec), ...
    channelName);

disp('Writing data to .DAT file ...');

if opts.Append && isfile(filePath)
    iAppend(filePath, data, hdr);
    return
end

h = spatialSlabIO('create', filePath, hdr);
cleanupObj = onCleanup(@() spatialSlabIO('close', h));
spatialSlabIO('write', h, 1:dimSizes(2), data);
spatialSlabIO('finalize', h);
clear cleanupObj
end

function [dimNames, dimSizes] = iResolveLayout(dimNamesIn, data, filePath)
%IRESOLVELAYOUT Validate DimNames and return them with the axis sizes of DATA.

if isempty(dimNamesIn)
    error('Umitoolbox:saveData:missingDimNames', ...
        ['Failed to save "%s": numeric data needs ''DimNames'', the axis layout ' ...
         'of the array, for example {''Y'',''X'',''T''} or {''Y'',''X'',''T'',''E''}.'], ...
        filePath);
end

dimNames = cellstr(string(dimNamesIn));
dimNames = dimNames(:).';
if ~iIsValidLayout(dimNames)
    error('Umitoolbox:saveData:invalidDimNames', ...
        ['Failed to save "%s": DimNames {%s} is not a valid .dat layout. It must ' ...
         'start with ''Y'',''X'', followed by distinct axes among ''T'',''E'',''F'' ' ...
         'in that order.'], filePath, strjoin(dimNames, ','));
end
if ndims(data) > numel(dimNames)
    error('Umitoolbox:saveData:invalidDimNames', ...
        ['Failed to save "%s": the array has %d dimensions but DimNames {%s} ' ...
         'names only %d axes.'], filePath, ndims(data), strjoin(dimNames, ','), ...
        numel(dimNames));
end
dimSizes = double(size(data, 1:numel(dimNames)));
end

function tf = iIsValidLayout(dimNames)
%IISVALIDLAYOUT 'Y','X', then distinct axes among 'T','E','F' in slot order.

tf = numel(dimNames) >= 2 && numel(dimNames) <= 5 && ...
    strcmp(dimNames{1}, 'Y') && strcmp(dimNames{2}, 'X');
if ~tf
    return
end
[isKnown, slot] = ismember(dimNames(3:end), {'T', 'E', 'F'});
tf = all(isKnown) && all(diff(slot) > 0);
end

function iAppend(filePath, data, hdr)
%IAPPEND Append DATA to an existing headered file along its last axis.

if ~isDatWithHeader(filePath)
    error('Umitoolbox:saveData:appendToHeaderless', ...
        ['Cannot append to "%s": the file has no header. Save the complete ' ...
         'array instead.'], filePath);
end

existing = readDatHeader(filePath);
fileInfo = dir(filePath);
expectedBytes = existing.dataOffset + existing.expectedDataBytes;

if ~existing.writeComplete
    iAppendMismatch(filePath, 'the previous write did not complete');
end
if ~strcmp(existing.dataClass, hdr.dataClass)
    iAppendMismatch(filePath, sprintf('the file stores %s values', existing.dataClass));
end
if ~isequal(existing.dimNames, hdr.dimNames)
    iAppendMismatch(filePath, sprintf('the file axes are {%s}, not {%s}', ...
        strjoin(existing.dimNames, ','), strjoin(hdr.dimNames, ',')));
end
nAxes = numel(existing.dimNames);
if nAxes < 3
    iAppendMismatch(filePath, sprintf(['the file axes {%s} have no axis after X ' ...
        'to grow'], strjoin(existing.dimNames, ',')));
end
for k = 1:nAxes - 1
    if existing.dimSizes(k) ~= hdr.dimSizes(k)
        iAppendMismatch(filePath, sprintf('axis %s is %d in the file but %d in the data', ...
            existing.dimNames{k}, existing.dimSizes(k), hdr.dimSizes(k)));
    end
end
if any(strcmp(existing.dimNames, 'T')) && single(hdr.frameRateHz) ~= single(existing.frameRateHz)
    iAppendMismatch(filePath, sprintf('the file frame rate is %g Hz, not %g Hz', ...
        existing.frameRateHz, hdr.frameRateHz));
end
if fileInfo.bytes ~= expectedBytes
    iAppendMismatch(filePath, 'the file size does not match its header');
end

grown = struct('dataClass', existing.dataClass, 'frameRateHz', existing.frameRateHz, ...
    'exposureMsec', existing.exposureMsec, 'channelName', existing.channelName, ...
    'dimNames', {existing.dimNames}, 'dimSizes', existing.dimSizes, ...
    'writeComplete', false);

fid = fopen(filePath, 'r+', 'ieee-le');
if fid == -1
    error('Umitoolbox:saveData:fileOpenFailed', ...
        'Could not open file for writing: "%s".', filePath);
end
cleanupObj = onCleanup(@() safeFclose(fid));

% Mark the file incomplete, append the data, then record the new size.
fwrite(fid, encodeDatHeader(grown), 'uint8');
fseek(fid, 0, 'eof');
nWritten = fwrite(fid, data, 'single');
if nWritten ~= numel(data)
    error('Umitoolbox:saveData:fileWriteFailed', ...
        'Failed to write all data to file "%s" (%d/%d elements written).', ...
        filePath, nWritten, numel(data));
end

grown.dimSizes(end) = grown.dimSizes(end) + hdr.dimSizes(end);
grown.writeComplete = true;
fseek(fid, 0, 'bof');
fwrite(fid, encodeDatHeader(grown), 'uint8');
clear cleanupObj
end

function iAppendMismatch(filePath, reason)
error('Umitoolbox:saveData:appendMismatch', ...
    'Cannot append to "%s": %s.', filePath, reason);
end

function [rate, exposure] = iResolveRateAndExposure(opts, filePath, hasT)
%IRESOLVERATEANDEXPOSURE FrameRateHz option, then Info; NaN without a T axis.

exposure = iValidExposure(opts.Info);

if ~hasT
    rate = NaN;
    return
end

if ~isempty(opts.FrameRateHz)
    rate = double(opts.FrameRateHz);
    return
end

rate = iValidRate(opts.Info, 'frameRateHz');
if isempty(rate)
    error('Umitoolbox:saveData:missingFrameRate', ...
        ['Failed to save "%s": the data has a T axis but no frame rate. Pass ' ...
         '''FrameRateHz'' or the source file''s ''Info'' (with frameRateHz).'], ...
        filePath);
end
end

function rate = iValidRate(S, fieldName)
rate = [];
if isstruct(S) && isscalar(S) && isfield(S, fieldName)
    value = S.(fieldName);
    if isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value) && value > 0
        rate = double(value);
    end
end
end

function exposure = iValidExposure(Info)
exposure = NaN;
if isstruct(Info) && isscalar(Info) && isfield(Info, 'exposureMsec')
    value = Info.exposureMsec;
    if isnumeric(value) && isscalar(value) && isreal(value)
        exposure = double(value);
    end
end
end

function save2umt(filePath, data)
%SAVE2UMT Save derived-data structure to a .umt MAT-file.

validateUMTStruct(data);

disp('Writing data to .UMT file ...');
save(filePath, '-struct', 'data', '-mat','-v7.3');
end
