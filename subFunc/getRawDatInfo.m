function [Ny, Nx, Nt, freqHz] = getRawDatInfo(SaveFolder, inFile)
%GETRAWDATINFO Return YXT dimensions (and frame rate) for a raw .dat file.
%
%   [Ny, Nx, Nt] = getRawDatInfo(SaveFolder, inFile)
%   [Ny, Nx, Nt, freqHz] = getRawDatInfo(SaveFolder, inFile)
%
% Prefer loadMetaData(...) so Nt follows the actual file size rather than a
% potentially stale AcqInfoStream.Length value.

if ~isfolder(SaveFolder)
    error('Umitoolbox:getRawDatInfo:MissingSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

meta = loadMetaData(inFile);

requiredFields = {'Height', 'Width', 'datLength'};
if nargout > 3
    requiredFields{end+1} = 'Freq';
end

if ~all(isfield(meta, requiredFields))
    error('Umitoolbox:getRawDatInfo:InvalidMetaData', ...
        'loadMetaData did not return %s for "%s".', strjoin(requiredFields, ', '), inFile);
end

Ny = double(meta.Height);
Nx = double(meta.Width);
Nt = double(meta.datLength);
if nargout > 3
    freqHz = double(meta.Freq);
end
end
