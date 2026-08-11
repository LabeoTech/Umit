function report = migrateLegacyOpticalLibrary(outputRoot, legacyRoot)
%MIGRATELEGACYOPTICALLIBRARY Convert legacy MAT optics to text/JSON resources.
%
%   report = migrateLegacyOpticalLibrary(outputRoot, legacyRoot)
%
%   Converts the production contents of SysSpect.mat, FilterSets.mat, and
%   CameraSpect.mat into the canonical UMIT Rig optical repertoire. Spectra
%   are written as two-column wavelength/response text files on 400:700 nm
%   and normalized by their maximum. Runtime optical lookup does not use the
%   legacy MAT files; this utility documents and tests the one-time migration.

if nargin < 1 || isempty(outputRoot)
    outputRoot = fullfile(fileparts(mfilename('fullpath')), 'resources', 'library');
end
if nargin < 2 || isempty(legacyRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    legacyRoot = fullfile(repoRoot, 'IOIAnalysis', 'SystemOpticalParams');
end

source = load(fullfile(legacyRoot, 'SysSpect.mat'));
legacyFilterSets = load(fullfile(legacyRoot, 'FilterSets.mat'));
legacyCameraAliases = load(fullfile(legacyRoot, 'CameraSpect.mat'));
wavelengthNm = (400:700)';

mapping = struct( ...
    'legacyID', {'Red', 'Green', 'Yellow', ...
        'FF01441511717', 'FF01496LP', 'FF01512630', 'FF0151444', ...
        'FF0154080', 'FF0155423', 'FF0156524', 'FF0161850', ...
        'FF0247230', 'FF02472SP', 'PF1024', 'PF1312', 'ThorQLux'}, ...
    'category', {'illumination', 'illumination', 'illumination', ...
        'filter', 'filter', 'filter', 'filter', 'filter', 'filter', ...
        'filter', 'filter', 'filter', 'filter', 'camera', 'camera', 'camera'}, ...
    'newID', {'LED_632nm', 'LED_521nm', 'LED_593nm', ...
        'FF01441511717', 'FF01496LP', 'FF01512630', 'FF0151444', ...
        'FF0154080', 'FF0155423', 'FF0156524', 'FF0161850', ...
        'FF0247230', 'FF02472SP', 'PF1024', 'PF1312', 'ThorQLux'});

writtenFiles = strings(0, 1);
for iMapping = 1:numel(mapping)
    item = mapping(iMapping);
    response = double(source.(item.legacyID)(:));
    if numel(response) ~= numel(wavelengthNm) || any(~isfinite(response)) || ...
            any(response < 0) || max(response) <= 0
        error('Umitoolbox:migrateLegacyOpticalLibrary:invalidLegacySpectrum', ...
            'Legacy spectrum "%s" is malformed.', item.legacyID);
    end
    response = response ./ max(response);
    destinationFolder = fullfile(outputRoot, 'spectra', item.category);
    if ~isfolder(destinationFolder)
        mkdir(destinationFolder);
    end
    destination = fullfile(destinationFolder, [item.newID '.txt']);
    metadata = struct( ...
        'id', item.newID, ...
        'legacyID', item.legacyID, ...
        'category', item.category, ...
        'displayName', item.legacyID, ...
        'manufacturer', '', ...
        'model', '', ...
        'peakWavelengthNm', wavelengthNm(find(response == max(response), 1)), ...
        'normalization', 'dividedByMaximum');
    localWriteSpectrum(destination, wavelengthNm, response, metadata);
    writtenFiles(end+1, 1) = string(destination); %#ok<AGROW>
end

filterSetFolder = fullfile(outputRoot, 'filterSets');
if ~isfolder(filterSetFolder)
    mkdir(filterSetFolder);
end
filterSetNames = fieldnames(legacyFilterSets);
for iSet = 1:numel(filterSetNames)
    legacyName = filterSetNames{iSet};
    legacyDefinition = legacyFilterSets.(legacyName);
    definition = struct( ...
        'id', legacyName, ...
        'displayName', legacyName, ...
        'excitationSpectrumID', localFilterReference(legacyDefinition.Excitation), ...
        'emissionSpectrumID', localFilterReference(legacyDefinition.Emission));
    destination = fullfile(filterSetFolder, [legacyName '.json']);
    localWriteJSON(destination, definition);
    writtenFiles(end+1, 1) = string(destination); %#ok<AGROW>
end

cameraAliases = struct();
cameraAliasNames = fieldnames(legacyCameraAliases);
for iAlias = 1:numel(cameraAliasNames)
    alias = cameraAliasNames{iAlias};
    cameraAliases.(alias) = legacyCameraAliases.(alias);
end
mappingFile = fullfile(outputRoot, 'legacyMappings.json');
localWriteJSON(mappingFile, struct( ...
    'illuminationAliases', struct('Amber', 'yellow'), ...
    'cameraAliases', cameraAliases, ...
    'spectra', mapping));
writtenFiles(end+1, 1) = string(mappingFile);

report = struct( ...
    'outputRoot', string(outputRoot), ...
    'spectrumCount', numel(mapping), ...
    'filterSetCount', numel(filterSetNames), ...
    'writtenFiles', writtenFiles);
end

function reference = localFilterReference(value)
if isempty(value) || strcmpi(char(string(value)), 'none')
    reference = '';
else
    reference = char(string(value));
end
end

function localWriteJSON(filePath, value)
fid = fopen(filePath, 'w');
if fid == -1
    error('Umitoolbox:migrateLegacyOpticalLibrary:fileWriteFailed', ...
        'Could not write file: %s', filePath);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(value, 'PrettyPrint', true));
end

function localWriteSpectrum(filePath, wavelengthNm, response, metadata)
fid = fopen(filePath, 'w');
if fid == -1
    error('Umitoolbox:migrateLegacyOpticalLibrary:fileWriteFailed', ...
        'Could not write file: %s', filePath);
end
cleanup = onCleanup(@() fclose(fid));
jsonLines = splitlines(string(jsonencode(metadata, 'PrettyPrint', true)));
for iLine = 1:numel(jsonLines)
    fprintf(fid, '# %s\n', char(jsonLines(iLine)));
end
fprintf(fid, '%d\t%.17g\n', [wavelengthNm(:), response(:)].');
end
