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
%   loadMetaData and opens the file for reading; it works for headered,
%   legacy sidecar, and AcqInfos-bound files whose first two axes are Y, X.
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

if ~isempty(varargin) && isstruct(varargin{1})
    switch lower(char(mode))
        case 'read'
            varargout{1} = iHandleRead(varargin{:});
        case 'close'
            iClose(varargin{:});
        otherwise
            error('Umitoolbox:spatialSlabIO:invalidInput', ...
                'Unknown mode "%s" for a spatialSlabIO handle. Use ''read'' or ''close''.', char(mode));
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
