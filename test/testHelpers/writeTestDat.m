function writeTestDat(filePath, data, frameRateHz, varargin)
%WRITETESTDAT Write a Y-X-T test .dat file (headered, or legacy sidecar).
%
%   writeTestDat(filePath, data, frameRateHz)
%   writeTestDat(filePath, data, frameRateHz, exposureMsec)
%   writeTestDat(..., 'Format', 'legacySidecar')
%   writeTestDat(..., 'DimNames', {'Y','X','T','E'})
%
%   .dat header Phase 5a. Shared by the tests that need .dat inputs.
%   DATA is a Y-by-X-by-T numeric array (T = 1 for a Y-by-X array); its
%   class is kept.
%
%   'Format':
%       'header' (default)  - headered file (spatialSlabIO create/write/
%                             finalize): axes Y,X,T, DATA's class, the given
%                             frame rate and exposure (NaN when omitted),
%                             channelName = the file's base name.
%                             'DimNames' (headered only; default Y,X,T)
%                             sets another layout, e.g. {'Y','X'} or
%                             {'Y','X','T','E'} (.dat header Phase 6b-1);
%                             the frame rate is NaN without a T axis.
%       'legacySidecar'     - headerless file plus a legacy sidecar
%                             <name>.mat (dim_names, datSize, datLength,
%                             Freq, Datatype), readable by loadMetaData
%                             without AcqInfos.mat. Used where a test needs
%                             headerless input.
%
%   Tests should not write headerless .dat files without a sidecar: those
%   are described only by AcqInfos.mat, which Phase 5b stops supporting.

p = inputParser;
addRequired(p, 'filePath', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'data', @(x) isnumeric(x) && ~isempty(x) && ndims(x) <= 5);
addRequired(p, 'frameRateHz', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addOptional(p, 'exposureMsec', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Format', 'header', @(x) any(strcmpi(x, {'header', 'legacySidecar'})));
addParameter(p, 'DimNames', {'Y', 'X', 'T'}, @(x) iscell(x) || isstring(x));
parse(p, filePath, data, frameRateHz, varargin{:});

filePath = char(p.Results.filePath);
[folder, base] = fileparts(filePath);
sizes = [size(data, 1), size(data, 2), size(data, 3)];
dataClass = class(data);

if strcmpi(p.Results.Format, 'legacySidecar')
    fid = fopen(filePath, 'w');
    assert(fid ~= -1, 'writeTestDat: cannot create "%s".', filePath);
    cleanupObj = onCleanup(@() fclose(fid));
    fwrite(fid, data, dataClass);
    clear cleanupObj
    sidecar = struct('dim_names', {{'Y', 'X', 'T'}}, 'datSize', sizes(1:2), ...
        'datLength', sizes(3), 'Freq', double(frameRateHz), 'Datatype', dataClass);
    save(fullfile(folder, [base '.mat']), '-struct', 'sidecar');
    return
end

dimNames = cellstr(string(p.Results.DimNames(:).'));
rate = double(frameRateHz);
if ~any(strcmp(dimNames, 'T'))
    rate = NaN;
end
info = struct('dataClass', dataClass, 'dimNames', {dimNames}, ...
    'dimSizes', size(data, 1:numel(dimNames)), ...
    'frameRateHz', rate, 'exposureMsec', double(p.Results.exposureMsec));
h = spatialSlabIO('create', filePath, datHeaderFromInfo(info, base));
cleanupObj = onCleanup(@() spatialSlabIO('close', h));
spatialSlabIO('write', h, 1:size(data, 2), data);
spatialSlabIO('finalize', h);
clear cleanupObj
end
