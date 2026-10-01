function varargout = spatialSlabIO(mode, varargin)
% SPATIALSLABIO  Stream spatial X-slabs from/to a .dat file.
%
% Handle-based reads (use these in new code):
%
%   h    = SPATIALSLABIO('open', filename)
%   h    = SPATIALSLABIO('open', filename, 'Info', Info)
%   slab = SPATIALSLABIO('read', h, xIdx)
%   slab = SPATIALSLABIO('read', h, xIdx, frameIdx)
%          SPATIALSLABIO('close', h)
%
%   'open' takes the file layout (data offset, class, shape) from
%   loadMetaData and opens the file for reading; it works for headered
%   and legacy sidecar files whose first two axes are Y, X.
%   The handle h is a struct with fields fid, filePath, Info, Ny, Nx,
%   trailingSizes (Info.dimSizes(3:end), 1 for a Y-X file), nFrames
%   (prod(trailingSizes)), dataClass, bytesPerValue, dataOffset.
%   With 'Info', open reuses an Info that loadMetaData already returned for
%   this same file (checked by path) instead of calling loadMetaData
%   again; use it when a file is reopened for many small reads.
%
%   'read' returns all Y for the columns xIdx, in the stored class:
%     - without frameIdx: all frames, size [Ny, numel(xIdx), trailingSizes];
%     - with frameIdx: linear indices into the trailing axes flattened in
%       memory order (1..nFrames), size [Ny, numel(xIdx), numel(frameIdx)].
%   Columns and frames are returned in the order given. Indices must be
%   unique integers in range.
%
%   'close' closes the file; closing twice is harmless.
%
% Handle-based writes (streamed outputs):
%
%   h = SPATIALSLABIO('create', filename, hdr)
%       SPATIALSLABIO('write', h, xIdx, slab)
%       SPATIALSLABIO('write', h, xIdx, slab, frameIdx)
%       SPATIALSLABIO('finalize', h)
%
%   'create' writes a version-1 header from HDR (the encodeDatHeader
%   fields: dataClass, frameRateHz, exposureMsec, channelName, dimNames,
%   dimSizes) with the write-complete flag at 0 and extends the file to its
%   final size without writing zeros. The handle has the same fields as an
%   'open' handle and also supports 'read'.
%
%   'write' mirrors 'read': SLAB is [Ny, numel(xIdx), trailingSizes] (all
%   frames) or [Ny, numel(xIdx), numel(frameIdx)]. Values are cast to the
%   file's class. One fwrite per run of contiguous columns and consecutive
%   frames, using the skip argument (Phase 3a benchmark, strategy W4).
%
%   'finalize' sets the write-complete flag and closes the file. 'close' on
%   a created handle leaves the file marked incomplete.
%
% Growable files (outputs whose length is only known at the end):
%
%   h = SPATIALSLABIO('create', filename, hdr, 'Growable', true)
%   h = SPATIALSLABIO('append', h, block)
%       SPATIALSLABIO('finalize', h)
%
%   A growable file has exactly three axes, Y and X first, and grows along
%   the last one. hdr.dimSizes(end) is ignored: the file starts with no
%   frames and no data bytes (its header size is 1 with the write-complete
%   flag clear; 0 would mean the axis is absent). 'append' writes BLOCK
%   ([Ny, Nx, nF], or [Ny, Nx] for one frame) after the last frame with one
%   fwrite, cast to the file's class, then rewrites the header with the
%   frames written so far and the flag still clear, so an interrupted
%   write leaves a readable, incomplete file. It returns the updated
%   handle (h.framesWritten). 'finalize' needs at least one frame and sets
%   the flag; 'close' leaves the file incomplete. 'read' and 'write' are
%   not available on growable handles, and 'append' only on them.
%
%   Reads use one fread per run of contiguous columns and consecutive
%   frames, with the fread skip argument jumping over the rest of each
%   frame (chosen from the Phase 3a slab I/O benchmark).
%
% Positional form (DEPRECATED for new code; kept unchanged for existing
% callers, including IOIAnalysis):
%
%   slab = SPATIALSLABIO('read', fid, Ny, Nx, Nt, xIdx, precision)
%   SPATIALSLABIO('write', fid, Ny, Nx, Nt, xIdx, precision, slab)
%
%   The caller opens the file; the binary file is assumed to store a
%   MATLAB-formatted array of size [Y,X,T] in column-major order starting
%   at byte 0.
%
%   Inputs:
%     mode      : 'read' or 'write'
%     fid       : file identifier from fopen
%     Ny, Nx, Nt: dimensions of on-disk data [Y,X,T]
%     xIdx      : vector of X indices to read/write
%     precision : numeric class (e.g. 'single', 'double')
%   Additional input (WRITE mode only):
%     slab      : array of size [Ny, numel(xIdx), Nt]
%   Output (READ mode only):
%     slab      : array of size [Ny, numel(xIdx), Nt]
%   Notes:
%     - If xIdx is contiguous, fast block I/O is used.
%     - Otherwise, a safe column-wise fallback is used.
%     - The file must already be preallocated in WRITE mode.
%
% See also: loadMetaData, loadData, mapDat, fread, fwrite, fseek

if (ischar(mode) || (isstring(mode) && isscalar(mode))) && strcmpi(mode, 'open')
    varargout{1} = iOpen(varargin{:});
    return
end

if (ischar(mode) || (isstring(mode) && isscalar(mode))) && strcmpi(mode, 'create')
    varargout{1} = iCreate(varargin{:});
    return
end

if ~isempty(varargin) && isstruct(varargin{1})
    switch lower(char(mode))
        case 'read'
            varargout{1} = iHandleRead(varargin{:});
        case 'write'
            iHandleWrite(varargin{:});
        case 'append'
            varargout{1} = iAppend(varargin{:});
        case 'finalize'
            iFinalize(varargin{:});
        case 'close'
            iClose(varargin{:});
        otherwise
            error('Umitoolbox:spatialSlabIO:invalidInput', ...
                ['Unknown mode "%s" for a spatialSlabIO handle. Use ''read'', ' ...
                 '''write'', ''append'', ''finalize'', or ''close''.'], char(mode));
    end
    return
end

if (ischar(mode) || (isstring(mode) && isscalar(mode))) && strcmpi(mode, 'close')
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        '''close'' needs a handle returned by spatialSlabIO(''open'', filename).');
end

% Deprecated positional form, unchanged.
if nargout > 0
    varargout{1} = iPositional(mode, varargin{:});
else
    iPositional(mode, varargin{:});
end
end

% =========================================================================
% Handle-based reads
% =========================================================================
function h = iOpen(filename, varargin)
if nargin < 1 || ~(ischar(filename) || (isstring(filename) && isscalar(filename)))
    error('Umitoolbox:spatialSlabIO:invalidInput', '''open'' needs a file name.');
end
filename = char(filename);
if isempty(varargin)
    Info = loadMetaData(filename);
elseif numel(varargin) == 2 && (ischar(varargin{1}) || isstring(varargin{1})) && ...
        strcmpi(varargin{1}, 'Info')
    Info = iCheckResolvedInfo(varargin{2}, filename);
else
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        '''open'' accepts only the name-value option ''Info''.');
end

if numel(Info.dimNames) < 2 || ~isequal(Info.dimNames(1:2), {'Y', 'X'})
    error('Umitoolbox:spatialSlabIO:unsupportedLayout', ...
        ['spatialSlabIO streams files whose first two axes are Y, X; "%s" has ' ...
         'axes {%s}. Use loadData instead.'], filename, strjoin(Info.dimNames, ','));
end

if strcmp(Info.format, 'header')
    machineFormat = 'ieee-le';
else
    machineFormat = 'n';
end
fid = fopen(Info.filePath, 'r', machineFormat);
if fid < 0
    error('Umitoolbox:spatialSlabIO:openFailed', 'Cannot open file for reading: %s', Info.filePath);
end

h = iMakeHandle(fid, Info, false, struct());
end

function h = iMakeHandle(fid, Info, writable, headerDescription)
% Common handle fields for 'open' and 'create'.
h = struct();
h.fid = fid;
h.filePath = fopen(fid);
h.Info = Info;
h.Ny = Info.dimSizes(1);
h.Nx = Info.dimSizes(2);
h.trailingSizes = Info.dimSizes(3:end);
if isempty(h.trailingSizes)
    h.trailingSizes = 1;
end
h.nFrames = prod(h.trailingSizes);
h.dataClass = Info.dataClass;
h.bytesPerValue = getByteSize(Info.dataClass);
h.dataOffset = Info.dataOffset;
h.writable = writable;
h.headerDescription = headerDescription;
h.growable = false;
h.framesWritten = 0;
end

% =========================================================================
% Handle-based writes
% =========================================================================
function h = iCreate(filename, hdr, varargin)
if nargin < 2 || ~(ischar(filename) || (isstring(filename) && isscalar(filename))) || ...
        ~isstruct(hdr) || ~isscalar(hdr)
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        '''create'' needs a file name and a header description struct.');
end
growable = false;
if ~isempty(varargin)
    if numel(varargin) == 2 && (ischar(varargin{1}) || isstring(varargin{1})) && ...
            strcmpi(varargin{1}, 'Growable') && islogical(varargin{2}) && isscalar(varargin{2})
        growable = varargin{2};
    else
        error('Umitoolbox:spatialSlabIO:invalidInput', ...
            '''create'' accepts only the name-value option ''Growable'' (logical scalar).');
    end
end
filename = char(filename);
if isempty(fileparts(filename))
    filename = fullfile(pwd, filename);
end
if growable
    h = iCreateGrowable(filename, hdr);
    return
end

% Validate and encode before touching the file: an invalid header creates nothing.
hdr.writeComplete = false;
headerBytes = encodeDatHeader(hdr);
described = decodeDatHeader(headerBytes);
if numel(described.dimNames) < 2 || ~isequal(described.dimNames(1:2), {'Y', 'X'})
    error('Umitoolbox:spatialSlabIO:unsupportedLayout', ...
        'spatialSlabIO creates files whose first two axes are Y, X; got {%s}.', ...
        strjoin(described.dimNames, ','));
end
totalBytes = described.dataOffset + described.expectedDataBytes;

fid = fopen(filename, 'w', 'ieee-le');
if fid < 0
    error('Umitoolbox:spatialSlabIO:openFailed', 'Cannot create file: %s', filename);
end
fwrite(fid, headerBytes, 'uint8');
fclose(fid);
iExtendFile(filename, totalBytes);

fid = fopen(filename, 'r+', 'ieee-le');
if fid < 0
    error('Umitoolbox:spatialSlabIO:openFailed', 'Cannot open file for writing: %s', filename);
end
% The new file is described from the header just written, so any file
% name works (for example a scratch "<name>.dat.tmp").
[Info, description] = iDescribeCreated(fid, described);
h = iMakeHandle(fid, Info, true, description);
end

function [Info, description] = iDescribeCreated(fid, described)
% Info (schema fields) and header description of a file being written.
Info = struct('filePath', fopen(fid), 'format', 'header', ...
    'dataOffset', described.dataOffset, 'dataClass', described.dataClass, ...
    'dimNames', {described.dimNames}, 'dimSizes', described.dimSizes, ...
    'frameRateHz', described.frameRateHz, 'exposureMsec', described.exposureMsec, ...
    'channelName', described.channelName, 'writeComplete', false);
description = struct('dataClass', described.dataClass, 'frameRateHz', described.frameRateHz, ...
    'exposureMsec', described.exposureMsec, 'channelName', described.channelName, ...
    'dimNames', {described.dimNames}, 'dimSizes', described.dimSizes);
end

function iExtendFile(filename, totalBytes)
% Grow the file to totalBytes without writing the data bytes.
info = dir(filename);
if info.bytes >= totalBytes
    return
end
try
    raf = java.io.RandomAccessFile(filename, 'rw');
    raf.setLength(totalBytes);
    raf.close();
catch
    % Without Java: zero-fill the remainder in chunks.
    fid = fopen(filename, 'a');
    cleanupObj = onCleanup(@() fclose(fid));
    remaining = totalBytes - info.bytes;
    chunk = 64 * 1024^2;
    while remaining > 0
        n = min(chunk, remaining);
        fwrite(fid, zeros(n, 1, 'uint8'), 'uint8');
        remaining = remaining - n;
    end
    clear cleanupObj
end
end

function iHandleWrite(h, xIdx, slab, frameIdx)
iAssertOpenHandle(h);
iAssertWritable(h);
iAssertNotGrowable(h, 'write');
if nargin < 3
    error('Umitoolbox:spatialSlabIO:invalidInput', '''write'' needs column indices and a slab.');
end
xIdx = iCheckIndices(xIdx, h.Nx, 'xIdx');
if nargin < 4
    frameIdx = 1:h.nFrames;
else
    frameIdx = iCheckIndices(frameIdx, h.nFrames, 'frameIdx');
end

Ny = h.Ny;
nX = numel(xIdx);
nF = numel(frameIdx);
if ~isnumeric(slab) && ~islogical(slab)
    error('Umitoolbox:spatialSlabIO:invalidInput', 'The slab must be numeric.');
end
if size(slab, 1) ~= Ny || size(slab, 2) ~= nX || numel(slab) ~= Ny * nX * nF
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        'The slab must have %d rows, %d columns, and %d frames.', Ny, nX, nF);
end
slab = reshape(cast(slab, h.dataClass), Ny, nX, nF);

frameBytes = h.Ny * h.Nx * h.bytesPerValue;
[colStart, colLen, colPos] = iRuns(xIdx);
[frmStart, frmLen, frmPos] = iRuns(frameIdx);
singleRun = isscalar(colStart) && isscalar(frmStart);

for c = 1:numel(colStart)
    nCols = colLen(c);
    skipBytes = (h.Nx - nCols) * Ny * h.bytesPerValue;
    blockPrecision = sprintf('%d*%s', Ny * nCols, h.dataClass);
    for f = 1:numel(frmStart)
        if singleRun
            block = slab;   % avoid copying the whole slab
        else
            block = slab(:, colPos{c}, frmPos{f});
        end
        offset = h.dataOffset + (frmStart(f) - 1) * frameBytes + ...
            (colStart(c) - 1) * Ny * h.bytesPerValue;
        if fseek(h.fid, offset, 'bof') ~= 0
            error('Umitoolbox:spatialSlabIO:writeFailed', 'Seek failed in %s.', h.filePath);
        end
        if skipBytes == 0
            written = fwrite(h.fid, block, h.dataClass);
        else
            % fwrite applies the skip before each block: write the first
            % frame's block, then the rest landing on the same columns.
            written = fwrite(h.fid, block(:, :, 1), h.dataClass);
            if frmLen(f) > 1
                written = written + fwrite(h.fid, block(:, :, 2:end), blockPrecision, skipBytes);
            end
        end
        if written ~= numel(block)
            error('Umitoolbox:spatialSlabIO:writeFailed', ...
                'Wrote %d of %d values to %s.', written, numel(block), h.filePath);
        end
    end
end
end

function iFinalize(h)
iAssertOpenHandle(h);
iAssertWritable(h);
description = h.headerDescription;
if isfield(h, 'growable') && h.growable
    if h.framesWritten < 1
        error('Umitoolbox:spatialSlabIO:invalidInput', ...
            'Cannot finalize "%s": no frames were appended.', h.filePath);
    end
    description.dimSizes(end) = h.framesWritten;
end
description.writeComplete = true;
headerBytes = encodeDatHeader(description);
if fseek(h.fid, 0, 'bof') ~= 0
    error('Umitoolbox:spatialSlabIO:writeFailed', 'Seek failed in %s.', h.filePath);
end
fwrite(h.fid, headerBytes, 'uint8');
fclose(h.fid);
end

function h = iCreateGrowable(filename, hdr)
% Growable file: header only, last-axis size 1 and write-complete clear.
if ~isfield(hdr, 'dimNames') || ~isfield(hdr, 'dimSizes') || ...
        numel(hdr.dimNames) ~= 3 || numel(hdr.dimSizes) ~= 3
    error('Umitoolbox:spatialSlabIO:unsupportedLayout', ...
        'A growable file must have exactly three axes (Y, X, and the axis it grows along).');
end
hdr.dimSizes(end) = 1;
hdr.writeComplete = false;
headerBytes = encodeDatHeader(hdr);
described = decodeDatHeader(headerBytes);
if ~isequal(described.dimNames(1:2), {'Y', 'X'})
    error('Umitoolbox:spatialSlabIO:unsupportedLayout', ...
        'spatialSlabIO creates files whose first two axes are Y, X; got {%s}.', ...
        strjoin(described.dimNames, ','));
end

fid = fopen(filename, 'w+', 'ieee-le');
if fid < 0
    error('Umitoolbox:spatialSlabIO:openFailed', 'Cannot create file: %s', filename);
end
fwrite(fid, headerBytes, 'uint8');

[Info, description] = iDescribeCreated(fid, described);
h = iMakeHandle(fid, Info, true, description);
h.growable = true;
h.framesWritten = 0;
end

function h = iAppend(h, block)
iAssertOpenHandle(h);
iAssertWritable(h);
if ~isfield(h, 'growable') || ~h.growable
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        '''append'' needs a handle from spatialSlabIO(''create'', ..., ''Growable'', true).');
end
if nargin < 2 || ~(isnumeric(block) || islogical(block)) || ndims(block) > 3 || ...
        size(block, 1) ~= h.Ny || size(block, 2) ~= h.Nx || isempty(block)
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        'The block must be numeric with %d rows and %d columns ([Ny, Nx, nF] or [Ny, Nx]).', ...
        h.Ny, h.Nx);
end
nF = size(block, 3);
frameBytes = h.Ny * h.Nx * h.bytesPerValue;
if fseek(h.fid, h.dataOffset + h.framesWritten * frameBytes, 'bof') ~= 0
    error('Umitoolbox:spatialSlabIO:writeFailed', 'Seek failed in %s.', h.filePath);
end
written = fwrite(h.fid, cast(block, h.dataClass), h.dataClass);
if written ~= numel(block)
    error('Umitoolbox:spatialSlabIO:writeFailed', ...
        'Wrote %d of %d values to %s.', written, numel(block), h.filePath);
end
h.framesWritten = h.framesWritten + nF;

% Record the frames written so far; the file stays marked incomplete.
description = h.headerDescription;
description.dimSizes(end) = h.framesWritten;
description.writeComplete = false;
if fseek(h.fid, 0, 'bof') ~= 0
    error('Umitoolbox:spatialSlabIO:writeFailed', 'Seek failed in %s.', h.filePath);
end
fwrite(h.fid, encodeDatHeader(description), 'uint8');
end

function iAssertNotGrowable(h, mode)
if isfield(h, 'growable') && h.growable
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        '''%s'' is not available on a growable handle; use ''append''.', mode);
end
end

function iAssertWritable(h)
if ~isfield(h, 'writable') || ~h.writable
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        'The handle for "%s" is read-only; write with a handle from spatialSlabIO(''create'', ...).', ...
        h.filePath);
end
end

function Info = iCheckResolvedInfo(Info, filename)
% Accept an Info already returned by loadMetaData for this same file.
schemaFields = {'filePath', 'format', 'dataOffset', 'dataClass', 'dimNames', 'dimSizes'};
if ~isstruct(Info) || ~isscalar(Info) || ~all(isfield(Info, schemaFields))
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        'Info must be the loadMetaData output for the file (missing .dat schema fields).');
end
if ~strcmp(iCanonicalPath(Info.filePath), iCanonicalPath(filename))
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        'Info describes "%s", not "%s".', Info.filePath, filename);
end
end

function p = iCanonicalPath(p)
p = char(p);
if isempty(fileparts(p))
    p = fullfile(pwd, p);
end
try
    p = char(java.io.File(p).getCanonicalPath());
catch
    % Without Java, compare the absolute paths as given.
end
end

function slab = iHandleRead(h, xIdx, frameIdx)
iAssertOpenHandle(h);
iAssertNotGrowable(h, 'read');
if nargin < 2
    error('Umitoolbox:spatialSlabIO:invalidInput', '''read'' needs column indices xIdx.');
end
xIdx = iCheckIndices(xIdx, h.Nx, 'xIdx');
allFrames = nargin < 3;
if allFrames
    frameIdx = 1:h.nFrames;
else
    frameIdx = iCheckIndices(frameIdx, h.nFrames, 'frameIdx');
end

Ny = h.Ny;
nX = numel(xIdx);
nF = numel(frameIdx);
slab = zeros(Ny, nX, nF, h.dataClass);

frameBytes = h.Ny * h.Nx * h.bytesPerValue;
[colStart, colLen, colPos] = iRuns(xIdx);
[frmStart, frmLen, frmPos] = iRuns(frameIdx);

for c = 1:numel(colStart)
    nCols = colLen(c);
    precision = sprintf('%d*%s=>%s', Ny * nCols, h.dataClass, h.dataClass);
    skipBytes = (h.Nx - nCols) * Ny * h.bytesPerValue;
    for f = 1:numel(frmStart)
        nFr = frmLen(f);
        offset = h.dataOffset + (frmStart(f) - 1) * frameBytes + ...
            (colStart(c) - 1) * Ny * h.bytesPerValue;
        if fseek(h.fid, offset, 'bof') ~= 0
            error('Umitoolbox:spatialSlabIO:readFailed', 'Seek failed in %s.', h.filePath);
        end
        values = fread(h.fid, Ny * nCols * nFr, precision, skipBytes);
        if numel(values) ~= Ny * nCols * nFr
            error('Umitoolbox:spatialSlabIO:readFailed', ...
                'Read %d of %d values from %s.', numel(values), Ny * nCols * nFr, h.filePath);
        end
        slab(:, colPos{c}, frmPos{f}) = reshape(values, Ny, nCols, nFr);
    end
end

if allFrames
    slab = reshape(slab, [Ny, nX, h.trailingSizes]);
end
end

function iClose(h)
if isscalar(h) && isfield(h, 'fid') && isfield(h, 'filePath') && iIsOpen(h)
    fclose(h.fid);
end
end

function iAssertOpenHandle(h)
requiredFields = {'fid', 'filePath', 'Ny', 'Nx', 'trailingSizes', 'nFrames', ...
    'dataClass', 'bytesPerValue', 'dataOffset'};
if ~isscalar(h) || ~all(isfield(h, requiredFields))
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        'The handle was not created by spatialSlabIO(''open'', filename).');
end
if ~iIsOpen(h)
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        'The handle for "%s" is closed.', h.filePath);
end
end

function tf = iIsOpen(h)
% The fid is still open on this handle's file (fids can be reused).
tf = isnumeric(h.fid) && isscalar(h.fid) && h.fid > 2 && ...
    strcmp(fopen(h.fid), h.filePath);
end

function idx = iCheckIndices(idx, n, name)
if ~isnumeric(idx) || isempty(idx) || ~isvector(idx) || ~isreal(idx) || ...
        any(~isfinite(idx)) || any(idx ~= round(idx)) || any(idx < 1) || any(idx > n)
    error('Umitoolbox:spatialSlabIO:invalidInput', ...
        '%s must be integers from 1 to %d.', name, n);
end
idx = double(idx(:).');
if numel(unique(idx)) ~= numel(idx)
    error('Umitoolbox:spatialSlabIO:invalidInput', '%s must not repeat indices.', name);
end
end

function [starts, lens, positions] = iRuns(idx)
% Split IDX (in the given order) into runs of consecutive increasing values.
breaks = [1, find(diff(idx) ~= 1) + 1, numel(idx) + 1];
nRuns = numel(breaks) - 1;
starts = zeros(1, nRuns);
lens = zeros(1, nRuns);
positions = cell(1, nRuns);
for r = 1:nRuns
    positions{r} = breaks(r):breaks(r + 1) - 1;
    starts(r) = idx(breaks(r));
    lens(r) = numel(positions{r});
end
end

% =========================================================================
% Deprecated positional form (unchanged behavior)
% =========================================================================
function slab = iPositional(mode, fid, Ny, Nx, Nt, xIdx, precision, slab)
    bytes = getByteSize(precision);
    nX = numel(xIdx);

    isRead  = strcmpi(mode, 'read');
    isWrite = strcmpi(mode, 'write');

    if ~(isRead || isWrite)
        error('spatialSlabIO:InvalidMode', ...
              'Mode must be ''read'' or ''write''.');
    end

    % ---------------------------------------------------------------------
    % Allocate / validate slab
    % ---------------------------------------------------------------------
    if isRead
        slab = zeros(Ny, nX, Nt, precision);
    else
        if nargin < 8
            error('spatialSlabIO:MissingInput', ...
                  'Slab input required in WRITE mode.');
        end

        if size(slab,1) ~= Ny || size(slab,2) ~= nX || size(slab,3) ~= Nt
            error('spatialSlabIO:SizeMismatch', ...
                  'Slab must be of size [Ny, numel(xIdx), Nt].');
        end
    end

    % ---------------------------------------------------------------------
    % Contiguous vs non-contiguous X indices
    % ---------------------------------------------------------------------
    isContiguous = all(diff(xIdx) == 1);

    if isContiguous
        % Fast path (block I/O)
        x0 = xIdx(1);

        for t = 1:Nt
            offset = ((t-1)*Nx*Ny + (x0-1)*Ny) * bytes;
            fseek(fid, offset, 'bof');

            if isRead
                slab(:,:,t) = fread(fid, [Ny, nX], precision);
            else
                fwrite(fid, slab(:,:,t), precision);
            end
        end

    else
        % Safe fallback (column-wise)
        for t = 1:Nt
            tOffset = (t-1) * Ny * Nx * bytes;

            for xi = 1:nX
                x = xIdx(xi);
                offset = tOffset + (x-1) * Ny * bytes;
                fseek(fid, offset, 'bof');

                if isRead
                    slab(:,xi,t) = fread(fid, Ny, precision);
                else
                    fwrite(fid, slab(:,xi,t), precision);
                end
            end
        end
    end
end
