function [outData, outFiles] = correctMotionArtifact(data, SaveFolder, varargin)
%CORRECTMOTIONARTIFACT Correct frame-wise movement artifacts in a YXT .dat file.
%
%   shifts = correctMotionArtifact(data, SaveFolder)
%   shifts = correctMotionArtifact(data, SaveFolder, Name, Value, ...)
%   [shifts, outFiles] = correctMotionArtifact(___)
%   info = correctMotionArtifact('pipelineInfo')
%
%   Estimates the frame-to-frame subpixel translation between every frame
%   of a reference Y-by-X-by-T .dat file and its first frame, using the
%   FFT-based DFTREGISTRATION algorithm (Guizar-Sicairos, Thurman & Fienup,
%   "Efficient subpixel image registration algorithms," Opt. Lett. 33,
%   156-158, 2008). The same per-frame shift is then applied to the
%   reference file and to every other .dat file found in SaveFolder, since
%   co-acquired channels/cameras share the same optical path and therefore
%   the same physical motion.
%
% Inputs
%   data       - Filename of the reference .dat file used to estimate the
%                per-frame shifts. If not found as given, it is resolved
%                relative to SaveFolder. Must be a continuous Y-X-T
%                single-precision file with at least 2 frames, and must be
%                one of the .dat files present in SaveFolder.
%   SaveFolder - Existing folder containing the reference file and every
%                other .dat file that should receive the same correction.
%                All .dat files in this folder must share the reference
%                file's Height, Width, and frame count.
%
% Name-Value Options
%   UpsamplingFactor - Positive integer upsampling factor passed to
%                      DFTREGISTRATION. 1 registers to the nearest whole
%                      pixel; higher values register to within 1/factor of
%                      a pixel. Default: 100.
%   Overwrite        - Logical scalar. When false (default), each
%                      corrected file is written next to the original as
%                      "<name>_MotionCorrected.dat" and no source .dat file
%                      is modified. When true, every .dat file in
%                      SaveFolder is destructively replaced in place.
%   SaveShifts       - Logical scalar. When true (default), the per-frame
%                      shift table and run metadata are saved to
%                      "MotionCorrectionShifts.mat" in SaveFolder.
%
% Output
%   outData  - Nt-by-2 double matrix of [rowShift, colShift] pixel shifts
%              estimated from the reference file (row 1 is always [0 0]).
%   outFiles - Cell array of full paths to every corrected .dat file that
%              was written, reference file first.
%
% Files Created or Modified
%   Overwrite=false (default): writes "<name>_MotionCorrected.dat" next to
%   every source .dat file in SaveFolder, plus "MotionCorrectionShifts.mat"
%   when SaveShifts is true. No source .dat file is modified.
%   Overwrite=true: every .dat file in SaveFolder is rewritten in place.
%
% Notes
%   - Pixels shifted in from outside the original frame are filled with 0.
%   - Raw .dat files are assumed to store continuous YXT data in single
%     precision, matching every other Analysis/ImageRegistration function.
%
% Example
%   shifts = correctMotionArtifact('acquisition_1.dat', saveFolder, ...
%       'UpsamplingFactor', 100);
%
% See also CREATEREGISTRATIONTFORM, APPLYREGISTRATIONTFORMONFOLDER.

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'data', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'UpsamplingFactor', 100, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == round(x));
addParameter(p, 'Overwrite', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SaveShifts', true, @(x) islogical(x) && isscalar(x));
parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
if ~isfolder(SaveFolder)
    error('Umitoolbox:correctMotionArtifact:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

usfac = double(p.Results.UpsamplingFactor);
overwrite = p.Results.Overwrite;
saveShifts = p.Results.SaveShifts;

refFile = char(string(p.Results.data));
if ~isfile(refFile)
    altPath = fullfile(SaveFolder, refFile);
    if isfile(altPath)
        refFile = altPath;
    else
        error('Umitoolbox:correctMotionArtifact:InputFileNotFound', ...
            'Reference file "%s" was not found.', data);
    end
end
[~, refName, refExt] = fileparts(refFile);
refFileName = [refName refExt];

datList = dir(fullfile(SaveFolder, '*.dat'));
datList = datList(~[datList.isdir]);
assert(~isempty(datList), ...
    'Umitoolbox:correctMotionArtifact:NoDatFiles', ...
    'No .dat files were found in "%s".', SaveFolder);

% Validate every target before estimating or applying any shift: this
% operation writes multiple files from a single set of per-frame shifts, so
% a size/layout mismatch discovered partway through would leave some
% channels corrected and others not.
[datPlan, refIdx] = iPreflightDatFiles(datList, SaveFolder, refFileName);

ny = datPlan(refIdx).ny;
nx = datPlan(refIdx).nx;
nt = datPlan(refIdx).nt;
if nt < 2
    error('Umitoolbox:correctMotionArtifact:TooFewFrames', ...
        'Reference file "%s" must contain at least 2 frames.', refFileName);
end

outData = iEstimateShifts(datPlan(refIdx).path, ny, nx, nt, usfac);

outFiles = cell(numel(datPlan), 1);
for iFile = 1:numel(datPlan)
    outFiles{iFile} = iApplyShiftsToFile(datPlan(iFile), outData, overwrite);
end

if saveShifts
    iSaveShiftsMat(SaveFolder, refFileName, usfac, outData, outFiles);
end

end

% =========================================================================
% Local helper: validate every .dat file against the reference layout
% =========================================================================
function [datPlan, refIdx] = iPreflightDatFiles(datList, SaveFolder, refFileName)
%IPREFLIGHTDATFILES Validate every .dat file against the reference layout.
%
% Every .dat file must share the reference file's Height/Width/frame count
% and use continuous Y-X-T layout: a per-frame shift vector is only valid
% when frame index t means the same acquisition instant in every file.

datPlan = struct('name', {}, 'path', {}, 'ny', {}, 'nx', {}, 'nt', {});

for iFile = 1:numel(datList)
    fileName = datList(iFile).name;
    datPath = fullfile(SaveFolder, fileName);
    md = loadMetaData(datPath);

    if ~all(isfield(md, {'Height','Width','datLength'}))
        error('Umitoolbox:correctMotionArtifact:InvalidMetadata', ...
            'Could not resolve Height/Width/datLength for "%s".', datPath);
    end

    dimNames = {'Y','X','T'};
    if isfield(md, 'dim_names') && ~isempty(md.dim_names)
        dimNames = cellstr(string(md.dim_names));
    end
    if ~isequal(dimNames(:).', {'Y','X','T'})
        error('Umitoolbox:correctMotionArtifact:UnsupportedLayout', ...
            ['File "%s" has dimensions {%s}. Motion correction only ' ...
             'supports continuous Y-X-T .dat files.'], ...
            fileName, strjoin(dimNames(:).', ','));
    end

    datPlan(end+1) = struct( ...
        'name', fileName, 'path', datPath, ...
        'ny', double(md.Height), 'nx', double(md.Width), ...
        'nt', double(md.datLength)); %#ok<AGROW>
end

refIdx = find(strcmp({datPlan.name}, refFileName), 1);
if isempty(refIdx)
    error('Umitoolbox:correctMotionArtifact:ReferenceNotInFolder', ...
        'Reference file "%s" must be one of the .dat files in "%s".', ...
        refFileName, SaveFolder);
end
refSize = [datPlan(refIdx).ny, datPlan(refIdx).nx, datPlan(refIdx).nt];

for iFile = 1:numel(datPlan)
    thisSize = [datPlan(iFile).ny, datPlan(iFile).nx, datPlan(iFile).nt];
    if ~isequal(thisSize, refSize)
        error('Umitoolbox:correctMotionArtifact:SizeMismatch', ...
            ['File "%s" is %dx%dx%d, but reference file "%s" is %dx%dx%d. ' ...
             'All .dat files in the folder must share the same frame size ' ...
             'and frame count to receive the same per-frame shift.'], ...
            datPlan(iFile).name, thisSize(1), thisSize(2), thisSize(3), ...
            refFileName, refSize(1), refSize(2), refSize(3));
    end
end

end

% =========================================================================
% Local helper: estimate per-frame shifts against frame 1
% =========================================================================
function shifts = iEstimateShifts(refPath, ny, nx, nt, usfac)
%IESTIMATESHIFTS Estimate per-frame DFTREGISTRATION shifts against frame 1.

frameBytes = ny * nx * getByteSize('single');

fid = fopen(refPath, 'r');
assert(fid ~= -1, 'Umitoolbox:correctMotionArtifact:FileOpenFailed', ...
    'Could not open "%s" for reading.', refPath);
c = onCleanup(@() safeFclose(fid));

frame1 = iReadFrame(fid, 0, ny, nx, frameBytes, refPath, 1);
refFFT = fft2(double(frame1));

shifts = zeros(nt, 2);
for t = 2:nt
    frame = iReadFrame(fid, t-1, ny, nx, frameBytes, refPath, t);
    regOutput = dftregistration(refFFT, fft2(double(frame)), usfac);
    shifts(t, :) = regOutput(3:4);
end

end

function frame = iReadFrame(fid, frameIdx0, ny, nx, frameBytes, filePath, frameNum)
%IREADFRAME Read one Y-by-X frame (0-based frameIdx0) as single precision.

fseek(fid, frameIdx0 * frameBytes, 'bof');
raw = fread(fid, ny * nx, '*single');
assert(numel(raw) == ny * nx, ...
    'Umitoolbox:correctMotionArtifact:InvalidDatFile', ...
    'Could not read frame %d from "%s".', frameNum, filePath);
frame = reshape(raw, ny, nx);

end

% =========================================================================
% Local helper: apply the estimated shifts to one .dat file
% =========================================================================
function outPath = iApplyShiftsToFile(planEntry, shifts, overwrite)
%IAPPLYSHIFTSTOFILE Apply the estimated per-frame shifts to one .dat file.

ny = planEntry.ny;
nx = planEntry.nx;
nt = planEntry.nt;
srcPath = planEntry.path;
frameBytes = ny * nx * getByteSize('single');

if overwrite
    destPath = srcPath;
else
    [srcFolder, stem, ext] = fileparts(srcPath);
    destPath = fullfile(srcFolder, [stem '_MotionCorrected' ext]);
end

[destFolder, destStem, destExt] = fileparts(destPath);
tmpPath = fullfile(destFolder, [destStem '_writing' destExt]);

fidIn = fopen(srcPath, 'r');
assert(fidIn ~= -1, 'Umitoolbox:correctMotionArtifact:FileOpenFailed', ...
    'Could not open "%s" for reading.', srcPath);
cIn = onCleanup(@() safeFclose(fidIn));

fidOut = fopen(tmpPath, 'w');
assert(fidOut ~= -1, 'Umitoolbox:correctMotionArtifact:FileOpenFailed', ...
    'Could not open temporary output file "%s".', tmpPath);
cOut = onCleanup(@() safeFclose(fidOut));

for t = 1:nt
    fseek(fidIn, (t-1) * frameBytes, 'bof');
    raw = fread(fidIn, ny * nx, '*single');
    assert(numel(raw) == ny * nx, ...
        'Umitoolbox:correctMotionArtifact:InvalidDatFile', ...
        'Could not read frame %d from "%s".', t, srcPath);
    frame = reshape(raw, ny, nx);

    rowShift = shifts(t, 1);
    colShift = shifts(t, 2);
    if rowShift ~= 0 || colShift ~= 0
        frame = imtranslate(frame, [colShift, rowShift], 'cubic', 'FillValues', 0);
    end

    fwrite(fidOut, frame, 'single');
end

clear cIn cOut % close both fids via safeFclose before the file move below

if overwrite
    iReplaceFileSafely(tmpPath, destPath);
else
    [moveOk, moveMsg] = movefile(tmpPath, destPath, 'f');
    assert(moveOk, 'Umitoolbox:correctMotionArtifact:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', tmpPath, destPath, moveMsg);
end

outPath = destPath;

end

% =========================================================================
% Local helper: atomic in-place replace (Overwrite=true only)
% =========================================================================
function iReplaceFileSafely(tmpPath, destPath)
%IREPLACEFILESAFELY Install tmpPath as destPath without a data-loss window.
%
% Mirrors applyRegistrationTformOnFolder's iReplaceFileSafely: move the
% existing file aside, install the replacement, then drop the backup, so a
% failed move cannot leave the folder without the data.

backupPath = [destPath '.bak'];
backupCreated = false;

if isfile(destPath)
    [ok, message] = movefile(destPath, backupPath, 'f');
    if ~ok
        error('Umitoolbox:correctMotionArtifact:BackupFailed', ...
            'Could not back up "%s" before replacing it: %s', destPath, message);
    end
    backupCreated = true;
end

[ok, message] = movefile(tmpPath, destPath, 'f');
if ~ok
    if backupCreated && isfile(backupPath)
        movefile(backupPath, destPath, 'f');
    end
    error('Umitoolbox:correctMotionArtifact:ReplaceFailed', ...
        'Could not install the corrected data as "%s": %s', destPath, message);
end

if backupCreated && isfile(backupPath)
    delete(backupPath);
end

end

% =========================================================================
% Local helper: persist shift table for provenance/QC
% =========================================================================
function iSaveShiftsMat(SaveFolder, refFileName, usfac, shifts, outFiles)
%ISAVESHIFTSMAT Persist the estimated shifts and run metadata for provenance.

MotionCorrection = struct( ...
    'referenceFile', refFileName, ...
    'upsamplingFactor', usfac, ...
    'shifts', shifts, ...
    'correctedFiles', {outFiles}, ...
    'appliedOn', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'appliedBy', mfilename);

save(fullfile(SaveFolder, 'MotionCorrectionShifts.mat'), 'MotionCorrection');

end

% =========================================================================
% Local pipeline info
% =========================================================================
function info = localPipelineInfo()
%LOCALPIPELINEINFO Return PipelineManager metadata for correctMotionArtifact.

info = PipelineManager.createPipelineInfo( ...
    'correctMotionArtifact', ...
    ['Estimate per-frame motion shifts from a reference .dat file and ' ...
     'apply them to every .dat file in the folder.']);

info = PipelineManager.addInput( ...
    info, ...
    'data', ...
    'ImageTimeSeries', ...
    'Reference Y-X-T .dat filename used to estimate per-frame shifts.', ...
    'kind', 'input', ...
    'position', 1, ...
    'callType', 'positional', ...
    'isData', true, ...
    'supportsFile', true, ...
    'dataMode', 'file');

info = PipelineManager.addInput( ...
    info, ...
    'SaveFolder', ...
    'SaveFolder', ...
    'Folder containing the reference file and every .dat file to correct.', ...
    'kind', 'input', ...
    'position', 2, ...
    'callType', 'positional', ...
    'isData', false);

info = PipelineManager.addInput( ...
    info, ...
    'UpsamplingFactor', ...
    'parameter', ...
    'DFTREGISTRATION upsampling factor (1 = whole-pixel, higher = subpixel).', ...
    'kind', 'parameter', ...
    'position', 3, ...
    'callType', 'namevalue', ...
    'default', 100, ...
    'allowed', [1 Inf], ...
    'dataType', 'double');

info = PipelineManager.addInput( ...
    info, ...
    'Overwrite', ...
    'parameter', ...
    ['If true, destructively rewrites every .dat file in place instead ' ...
     'of writing "_MotionCorrected" copies.'], ...
    'kind', 'parameter', ...
    'position', 4, ...
    'callType', 'namevalue', ...
    'default', false, ...
    'allowed', [true false], ...
    'dataType', 'logical');

info = PipelineManager.addInput( ...
    info, ...
    'SaveShifts', ...
    'parameter', ...
    'If true, saves the per-frame shift table to MotionCorrectionShifts.mat.', ...
    'kind', 'parameter', ...
    'position', 5, ...
    'callType', 'namevalue', ...
    'default', true, ...
    'allowed', [true false], ...
    'dataType', 'logical');

info = PipelineManager.addOutput( ...
    info, ...
    'outData', ...
    'MotionShifts', ...
    'data', ...
    'Nt-by-2 [rowShift, colShift] pixel shifts estimated from the reference file.', ...
    'MotionCorrectionShifts.mat', ...
    1, ...
    'isData', true, ...
    'saveFileName', 'MotionCorrectionShifts.mat');

info = PipelineManager.addOutput( ...
    info, ...
    'correctedDatFiles', ...
    {'ImageTimeSeries','ProcessedData'}, ...
    'file', ...
    ['Every .dat file in SaveFolder, corrected and written as a new file ' ...
     '(or in place when Overwrite is true).'], ...
    '*.dat', ...
    2, ...
    'isData', false, ...
    'returnsValue', false);

info.notes = { ...
    ['Uses the third-party DFTREGISTRATION algorithm (Guizar-Sicairos, ' ...
     'Thurman & Fienup, Opt. Lett. 33, 156-158, 2008), bundled as a ' ...
     'private helper.']; ...
    'Overwrite=false (default) never modifies a source .dat file.'};

end
