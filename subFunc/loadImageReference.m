function [ImageReference, referenceFileAbs] = loadImageReference(referenceFile)
%LOADIMAGEREFERENCE Load and validate an ImageReference file.
%
%   ImageReference = loadImageReference(referenceFile)
%   [ImageReference, referenceFileAbs] = loadImageReference(referenceFile)
%
%   referenceFile can be:
%       - an absolute path to ImageReference_*.mat
%       - a path relative to the UMIT root folder

arguments
    referenceFile (1,1) string
end

referenceFile = strtrim(referenceFile);

if referenceFile == ""
    error('loadImageReference:MissingFile', ...
        'referenceFile must be provided.');
end

if isfile(referenceFile)
    referenceFileAbs = char(referenceFile);
else
    referenceFileAbs = fullfile(getUmitFolder(), char(referenceFile));
end

if ~isfile(referenceFileAbs)
    error('loadImageReference:FileNotFound', ...
        'ImageReference file not found: %s', referenceFileAbs);
end

S = load(referenceFileAbs, 'ImageReference');

if ~isfield(S, 'ImageReference')
    error('loadImageReference:MissingVariable', ...
        'File does not contain variable "ImageReference": %s', referenceFileAbs);
end

ImageReference = S.ImageReference;
validateImageReference(ImageReference);

end