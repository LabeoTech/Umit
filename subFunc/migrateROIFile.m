function ROIFile = migrateROIFile(ROIFile)
%MIGRATEROIFILE Upgrade ROIFile to the current schema version.
%
%   ROIFile = migrateROIFile(ROIFile)
%
%   Version 1 is the first versioned ROI schema. Pre-version legacy ROI
%   structures are intentionally not supported yet.

if ~isstruct(ROIFile) || ~isscalar(ROIFile)
    error('migrateROIFile:InvalidROIFile', ...
        'ROIFile must be a scalar structure.');
end

if ~isfield(ROIFile, 'version')
    error('migrateROIFile:MissingVersion', ...
        'ROI file does not contain a version field. Retrocompatibility will be added later.');
end

switch ROIFile.version
    case 1
        % Current schema. No migration needed.
        return

    otherwise
        error('migrateROIFile:UnsupportedVersion', ...
            'Unsupported ROI file version: %g.', ROIFile.version);
end

end
