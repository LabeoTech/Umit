function outFile = saveData(filename, data, varargin)
%SAVEDATA Save raw image-series data or derived UMIT data to disk.
%
%   saveData(filename, data)
%   saveData(filename, data, AcqInfoStream)
%   saveData(..., 'Info', Info)
%   saveData(..., 'FrameRateHz', rate)
%   saveData(..., 'ChannelName', name)
%   saveData(..., 'Append', true)
%
%   Inputs:
%       filename        - Full path of the file to be saved. The extension
%                         is forced based on the input data type.
%       data            - Either:
%                         1) non-empty single numeric 3D array (Y-X-T),
%                            saved as a headered .dat file
%                         2) valid derived-data structure, saved as .umt
%
%   Optional inputs for numeric data only:
%       AcqInfoStream   - Acquisition info structure. Saved as AcqInfos.mat
%                         when that file does not exist yet, and used as
%                         the frame-rate fallback (see Notes).
%
%   Name-Value options (numeric data only):
%       Info            - .dat Info struct of the source file (as returned
%                         by loadMetaData). Its frameRateHz and exposureMsec
%                         fields, when present, describe the output.
%       FrameRateHz     - Positive finite frame rate of the output.
%                         Overrides Info.frameRateHz.
%       ChannelName     - Header channelName. Default: the output file's
%                         base name. (PipelineManager passes the save
%                         target's name for temporary files that may be
%                         promoted to it.)
%       Append          - Logical scalar. If true, append the frames to an
%                         existing headered .dat file along its last axis.
%                         Default: false
%
%   Notes:
%       - Numeric arrays are written with a version-1 header (see
%         docs/dev/dat-header-spec.md): axes {'Y','X','T'}, class single,
%         channelName = the output file's base name (non-ASCII characters
%         replaced by '_', truncated to the header field).
%       - Frame rate: FrameRateHz, else Info.frameRateHz, else the rate
%         "AcqInfos.mat" in the target folder (or the AcqInfoStream
%         argument) gives this file: its ImportedChannels entry matched by
%         file name, then by length (resolveDatTimeline, as headerless
%         loading does), else its top-level FrameRateHz. Error if none gives
%         a finite rate > 0. The AcqInfos fallback is transitional.
%       - Exposure: Info.exposureMsec; else, when the rate came from an
%         ImportedChannels entry, that entry's ExposureMsec; else NaN.
%       - An existing file is overwritten, headered or not.
%       - Append to a missing file creates it. Append to a headered file
%         requires the same class, the same Y and X sizes, a Y-X-T layout,
%         a complete earlier write, and the same frame rate. Append to a
%         headerless file is an error and leaves the file unchanged.
%       - Structured derived data are saved as MAT-files with ".umt".

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'saveData';

addRequired(p, 'filename', @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(p, 'data', @(x) (isnumeric(x) && isa(x, 'single') && ~isempty(x)) || isstruct(x));
addOptional(p, 'AcqInfoStream', struct.empty(0,1), @isstruct);
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
%SAVE2DAT Save a numeric Y-X-T array to a headered .dat file.

if ~(isnumeric(data) && isa(data, 'single') && ~isempty(data))
    error('Umitoolbox:saveData:invalidInput', ...
        'Numeric data must be a non-empty single array.');
end

if ndims(data) ~= 3
    error('Umitoolbox:saveData:invalidInput', ...
        'Numeric data must be a 3D single array.');
end

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

acqInfoFile = fullfile(saveFolder, 'AcqInfos.mat');
if ~isfile(acqInfoFile) && ~isempty(opts.AcqInfoStream)
    AcqInfoStream = opts.AcqInfoStream;
    save(acqInfoFile, 'AcqInfoStream');
end

[frameRateHz, exposureMsec] = iResolveRateAndExposure(opts, acqInfoFile, filePath, size(data, 3));

channelName = fileName;
if strlength(string(opts.ChannelName)) > 0
    channelName = char(opts.ChannelName);
end
hdr = datHeaderFromInfo(struct('dataClass', 'single', 'dimNames', {{'Y', 'X', 'T'}}, ...
    'dimSizes', size(data), 'frameRateHz', frameRateHz, 'exposureMsec', exposureMsec), ...
    channelName);

disp('Writing data to .DAT file ...');

if opts.Append && isfile(filePath)
    iAppend(filePath, data, hdr);
    return
end

h = spatialSlabIO('create', filePath, hdr);
cleanupObj = onCleanup(@() spatialSlabIO('close', h));
spatialSlabIO('write', h, 1:size(data, 2), data);
spatialSlabIO('finalize', h);
clear cleanupObj
end

function iAppend(filePath, data, hdr)
%IAPPEND Append frames to an existing headered file along its last axis.

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
    iAppendMismatch(filePath, sprintf('the file axes are {%s}', strjoin(existing.dimNames, ',')));
end
if ~isequal(existing.dimSizes(1:2), hdr.dimSizes(1:2))
    iAppendMismatch(filePath, sprintf('the file frames are %d x %d', ...
        existing.dimSizes(1), existing.dimSizes(2)));
end
if single(hdr.frameRateHz) ~= single(existing.frameRateHz)
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

% Mark the file incomplete, append the frames, then record the new size.
fwrite(fid, encodeDatHeader(grown), 'uint8');
fseek(fid, 0, 'eof');
nWritten = fwrite(fid, data, 'single');
if nWritten ~= numel(data)
    error('Umitoolbox:saveData:fileWriteFailed', ...
        'Failed to write all data to file "%s" (%d/%d elements written).', ...
        filePath, nWritten, numel(data));
end

grown.dimSizes(end) = grown.dimSizes(end) + size(data, 3);
grown.writeComplete = true;
fseek(fid, 0, 'bof');
fwrite(fid, encodeDatHeader(grown), 'uint8');
clear cleanupObj
end

function iAppendMismatch(filePath, reason)
error('Umitoolbox:saveData:appendMismatch', ...
    'Cannot append to "%s": %s.', filePath, reason);
end

function [rate, exposure] = iResolveRateAndExposure(opts, acqInfoFile, filePath, nFrames)
%IRESOLVERATEANDEXPOSURE FrameRateHz option, then Info, then AcqInfos.mat.

exposure = iValidExposure(opts.Info);

if ~isempty(opts.FrameRateHz)
    rate = double(opts.FrameRateHz);
    return
end

rate = iValidRate(opts.Info, 'frameRateHz');
if ~isempty(rate)
    return
end

if isfile(acqInfoFile)
    S = load(acqInfoFile);
    if isfield(S, 'AcqInfoStream') && isstruct(S.AcqInfoStream) && isscalar(S.AcqInfoStream)
        [rate, channelExposure] = iRateFromAcqInfos(S.AcqInfoStream, filePath, nFrames);
        if isnan(exposure)
            exposure = channelExposure;
        end
    end
end
if isempty(rate)
    error('Umitoolbox:saveData:missingFrameRate', ...
        ['Failed to save "%s": no frame rate. Pass ''FrameRateHz'' or the ' ...
         'source file''s ''Info'', or provide AcqInfos.mat with FrameRateHz.'], ...
        filePath);
end
end

function [rate, exposure] = iRateFromAcqInfos(acq, filePath, nFrames)
%IRATEFROMACQINFOS Rate (and channel exposure) AcqInfos.mat gives this file.

rate = [];
exposure = NaN;
try
    timeline = resolveDatTimeline(nFrames, acq, 'DatFile', filePath, 'ThrowError', false);
catch
    timeline = struct('IsValid', false);
end

if timeline.IsValid
    rate = iValidRate(timeline, 'FrameRateHz');
    if ~isempty(rate) && startsWith(timeline.SourceType, 'ImportedChannels') && ...
            isfield(acq, 'ImportedChannels') && isfield(acq.ImportedChannels, 'ExposureMsec') && ...
            timeline.SourceIndex >= 1 && timeline.SourceIndex <= numel(acq.ImportedChannels)
        exposure = iValidExposure(struct('exposureMsec', ...
            acq.ImportedChannels(timeline.SourceIndex).ExposureMsec));
    end
end
if isempty(rate)
    rate = iValidRate(acq, 'FrameRateHz');
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
