function source = normalizeReferenceSourceInput(source)
%NORMALIZEREFERENCESOURCEINPUT Reference source struct from a struct or a path.
%
%   source = normalizeReferenceSourceInput(source)
%
%   A struct is returned unchanged (resolveImageReferenceSource checks it).
%   A text path becomes a 'umitDat' source for a .dat file and a 'file'
%   source otherwise, so pipelines can pass a plain path parameter.
%
%   See also resolveImageReferenceSource.

if isstruct(source)
    return
end
if ~(ischar(source) || (isstring(source) && isscalar(source))) || strlength(string(source)) == 0
    error('Umitoolbox:ImageRegistration:invalidSource', ...
        'A reference source is a source struct or a file path.');
end
path = char(string(source));
[~, ~, ext] = fileparts(path);
if strcmpi(ext, '.dat')
    source = struct('kind', 'umitDat', 'datPath', path);
else
    source = struct('kind', 'file', 'path', path);
end
end
