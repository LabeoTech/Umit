function [outData, outFiles] = correctMotionArtifact(data, SaveFolder, varargin)
%CORRECTMOTIONARTIFACT Correct frame-wise motion artifacts in a YXT .dat file.
%
%   [outData, outFiles] = correctMotionArtifact(data, SaveFolder)
%   [outData, outFiles] = correctMotionArtifact(data, SaveFolder, Name, Value, ...)
%   info = correctMotionArtifact('pipelineInfo')
%
%   Estimates the transform that registers every frame of a Y-by-X-by-T
%   .dat file to its own first frame, and warps each frame with it using
%   cubic interpolation. Frames are high-pass filtered before registration
%   so that uneven illumination does not bias the estimate.
%
%   Each file is corrected independently, against its own first frame. The
%   transforms are not shared between files: with strobed illumination the
%   channels are sampled at different instants, so motion can differ between
%   them. To correct several channels, run the function once per file.
%
%   Two transform types are available (TransformType):
%     'translation' - row/column shift estimated with the FFT-based
%                     DFTREGISTRATION algorithm (Guizar-Sicairos, Thurman &
%                     Fienup, "Efficient subpixel image registration
%                     algorithms," Opt. Lett. 33, 156-158, 2008). Fast.
%     'similarity'  - rotation, isotropic scale and translation estimated by
%                     phase correlation (IMREGCORR). About 5x slower.
%
%   When finished, a summary of the correlation of every frame with the
%   first frame, before and after correction, is printed and plotted.
%
% Inputs
%   data       - Filename of the .dat file to correct. If not found as
%                given, it is resolved relative to SaveFolder. Must be a
%                continuous Y-X-T single-precision file with at least 2
%                frames, located in SaveFolder.
%   SaveFolder - Existing folder containing the file. Corrected files and
%                run provenance are written here.
%
% Name-Value Options
%   TransformType    - 'translation' (default) or 'similarity'.
%   UpsamplingFactor - Positive integer upsampling factor of the
%                      'translation' estimate. 1 registers to the nearest
%                      whole pixel; higher values register to within
%                      1/factor of a pixel. Ignored for 'similarity'.
%                      Default: 100.
%   Overwrite        - Logical scalar. When false (default), the corrected
%                      file is written next to the original as
%                      "<name>_MotionCorrected.dat" and the source file is
%                      not modified. When true, the file is destructively
%                      replaced in place.
%   SaveShifts       - Logical scalar. When true (default), the per-frame
%                      transform parameters, correlations and run metadata
%                      are saved to "<name>_MotionCorrection.mat", and the
%                      QC figure to "<name>_MotionCorrectionQC.png", in
%                      SaveFolder.
%   ShowPlot         - Logical scalar. When true (default), shows the QC
%                      figure.
%
% Output
%   outData  - Nt-by-4 double matrix [tx, ty, rotationDeg, scale] of the
%              transform of each frame onto the first frame, in pixels and
%              degrees (tx is the column shift, ty the row shift; rotation 0
%              and scale 1 for 'translation'). Row 1 is always [0 0 0 1].
%   outFiles - Cell array with the full path of the corrected .dat file.
%
% Files Created or Modified
%   Overwrite=false (default): writes "<name>_MotionCorrected.dat", plus the
%   provenance .mat and QC figure when SaveShifts is true. The source file
%   is not modified.
%   Overwrite=true: the .dat file is rewritten in place.
%
% Notes
%   - Pixels warped in from outside the frame are filled with 0. NaN pixels
%     in the source are also set to 0 in the corrected file.
%   - Cubic interpolation applies sub-pixel shifts without rounding, which
%     avoids the whole-pixel flips that nearest-neighbor interpolation
%     introduces between consecutive frames. It slightly smooths the image,
%     so corrected pixel values are interpolated, not original.
%   - If the transform of a frame cannot be estimated, a warning is issued
%     and that frame is left uncorrected (identity transform).
%   - The reported correlation is the Pearson correlation between the
%     high-pass filtered first frame and each filtered frame, computed on
%     the image excluding a 5% border so that zero-filled edges do not count
%     against the correction. The summary counts frames whose correlation
%     dropped by more than 0.01.
%
% Example
%   [shifts, outFiles] = correctMotionArtifact('green.dat', saveFolder, ...
%       'TransformType', 'similarity');
%
% See also CREATEREGISTRATIONFORM, APPLYREGISTRATIONFORMONFOLDER.

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'data', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'TransformType', 'translation', ...
    @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    ismember(lower(char(string(x))), {'translation','similarity'}));
addParameter(p, 'UpsamplingFactor', 100, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == round(x));
addParameter(p, 'Overwrite', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SaveShifts', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowPlot', true, @(x) islogical(x) && isscalar(x));
parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
if ~isfolder(SaveFolder)
    error('Umitoolbox:correctMotionArtifact:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

transformType = lower(char(string(p.Results.TransformType)));
usfac = double(p.Results.UpsamplingFactor);
overwrite = p.Results.Overwrite;
saveShifts = p.Results.SaveShifts;
showPlot = p.Results.ShowPlot;

% Show warnings as a single line, without the call stack.
warnState = warning('off', 'backtrace');
restoreWarnings = onCleanup(@() warning(warnState));

inFile = char(string(p.Results.data));
if ~isfile(inFile)
    altPath = fullfile(SaveFolder, inFile);
    if isfile(altPath)
        inFile = altPath;
    else
        error('Umitoolbox:correctMotionArtifact:InputFileNotFound', ...
            'Input file "%s" was not found.', data);
    end
end
[~, inStem, inExt] = fileparts(inFile);
inFileName = [inStem inExt];

% The corrected file and provenance are written beside the source, so the
% source must live in SaveFolder.
iAssertInFolder(inFile, SaveFolder, inFileName);

entry = iResolveDatFile(fullfile(SaveFolder, inFileName), inFileName);
if entry.nt < 2
    error('Umitoolbox:correctMotionArtifact:TooFewFrames', ...
        'File "%s" must contain at least 2 frames.', inFileName);
end

fprintf('Estimating %s motion in %s...\n', transformType, inFileName);
[tforms, outData, failedFrames, corrBefore] = iEstimateTransforms(entry, transformType, usfac);

fprintf('Correcting %s...\n', inFileName);
[outPath, corrAfter] = iApplyTransformsToFile(entry, tforms, overwrite);
outFiles = {outPath};

iPrintSummary(entry, outData, corrBefore, corrAfter, failedFrames);

qcFile = '';
if saveShifts
    qcFile = fullfile(SaveFolder, [inStem '_MotionCorrectionQC.png']);
end
if saveShifts || showPlot
    iPlotQC(entry, transformType, outData, corrBefore, corrAfter, showPlot, qcFile);
end

if saveShifts
    iSaveProvenance(SaveFolder, entry, transformType, usfac, outData, ...
        corrBefore, corrAfter, failedFrames, outPath);
end

end

% =========================================================================
% Local helper: validate the input file
% =========================================================================
function iAssertInFolder(inFile, SaveFolder, inFileName)
%IASSERTINFOLDER Require that inFile lives directly in SaveFolder.

fileDir = dir(inFile);
folderDir = dir(fullfile(SaveFolder, '.'));
sameFolder = ~isempty(fileDir) && ~isempty(folderDir) && ...
    strcmpi(fileDir(1).folder, folderDir(1).folder);
if ~sameFolder
    error('Umitoolbox:correctMotionArtifact:InputNotInFolder', ...
        'Input file "%s" must be a .dat file located in "%s".', inFileName, SaveFolder);
end

end

function entry = iResolveDatFile(datPath, fileName)
%IRESOLVEDATFILE Read and validate the layout of the .dat file.
%
% Corrected frames are indexed by time, so the file must be a continuous
% Y-X-T single-precision movie.

md = loadMetaData(datPath);

if ~all(isfield(md, {'dimNames','dimSizes','dataClass'}))
    error('Umitoolbox:correctMotionArtifact:InvalidMetadata', ...
        'Could not resolve the axes, sizes, and class of "%s".', datPath);
end

dimNames = cellstr(string(md.dimNames));
if ~isequal(dimNames(:).', {'Y','X','T'})
    error('Umitoolbox:correctMotionArtifact:UnsupportedLayout', ...
        ['File "%s" has dimensions {%s}. Motion correction only ' ...
         'supports continuous Y-X-T .dat files.'], ...
        fileName, strjoin(dimNames(:).', ','));
end

if ~strcmp(md.dataClass, 'single')
    error('Umitoolbox:correctMotionArtifact:unsupportedDataClass', ...
        'Motion correction supports single-precision .dat files; "%s" stores %s.', ...
        fileName, md.dataClass);
end

[~, stem] = fileparts(fileName);
entry = struct('name', fileName, 'stem', stem, 'path', datPath, ...
    'ny', datAxisSize(md, 'Y'), 'nx', datAxisSize(md, 'X'), ...
    'nt', datAxisSize(md, 'T'), 'info', md);

end

% =========================================================================
% Local helpers: estimate the transform of every frame onto frame 1
% =========================================================================
function [tforms, params, failedFrames, corrBefore] = iEstimateTransforms(entry, transformType, usfac)
%IESTIMATETRANSFORMS Estimate the transform of each frame onto frame 1.
%
% Also returns the correlation of each (uncorrected) frame with frame 1.

slabIn = spatialSlabIO('open', entry.path, 'Info', entry.info);
c = onCleanup(@() spatialSlabIO('close', slabIn));

fixed = iPrepareForRegistration(iReadFrame(slabIn, entry, 1));
if strcmp(transformType, 'translation')
    fixedFFT = fft2(fixed);
end
tforms = cell(entry.nt, 1);
params = zeros(entry.nt, 4);
params(1, :) = [0 0 0 1];
tforms{1} = iMakeTform(eye(3));
corrBefore = ones(entry.nt, 1);
failedFrames = [];
for t = 2:entry.nt
    moving = iPrepareForRegistration(iReadFrame(slabIn, entry, t));
    corrBefore(t) = iFrameCorrelation(fixed, moving);
    try
        switch transformType
            case 'translation'
                reg = dftregistration(fixedFFT, fft2(moving), usfac);
                tx = reg(4);   % column shift
                ty = reg(3);   % row shift
                tforms{t} = iMakeTform([1 0 tx; 0 1 ty; 0 0 1]);
                params(t, :) = [tx ty 0 1];
            case 'similarity'
                tforms{t} = imregcorr(moving, fixed, 'similarity');
                params(t, :) = iTformParams(tforms{t});
        end
    catch ME
        warning('Umitoolbox:correctMotionArtifact:RegistrationFailed', ...
            'Frame %d of %s was left uncorrected: %s', t, entry.name, ME.message);
        tforms{t} = iMakeTform(eye(3));
        params(t, :) = [0 0 0 1];
        failedFrames(end+1) = t; %#ok<AGROW>
    end
    if mod(t, 500) == 0
        fprintf('  %s: %d/%d frames\n', entry.name, t, entry.nt);
    end
end

end

function frame = iReadFrame(slabIn, entry, t)
%IREADFRAME Read frame t (1-based) as a Y-by-X single array.

frame = reshape(spatialSlabIO('read', slabIn, 1:entry.nx, t), entry.ny, entry.nx);

end

function frame = iPrepareForRegistration(frame)
%IPREPAREFORREGISTRATION High-pass filter a frame before registration.
%
% Removes the slowly varying illumination (vignetting) that otherwise
% dominates the correlation and pulls the estimate towards no motion. The
% filter ignores NaN pixels (normalized convolution) and sets them to 0, so
% a masked border does not create an artificial edge.

frame = double(frame);
valid = ~isnan(frame);
frame(~valid) = 0;
w = double(valid);
sigSmooth = 0.5;
sigBackground = 0.05 * max(size(frame));
smoothed = imgaussfilt(frame, sigSmooth) ./ max(imgaussfilt(w, sigSmooth), eps);
background = imgaussfilt(frame, sigBackground) ./ max(imgaussfilt(w, sigBackground), eps);
frame = (smoothed - background) .* w;
frame = frame - mean(frame, 'all');

end

function r = iFrameCorrelation(a, b)
%IFRAMECORRELATION Pearson correlation of two prepared frames, excluding a 5% border.
%
% The border keeps the zero-filled edges created by a shift from counting
% against the correction.

mY = round(0.05 * size(a, 1));
mX = round(0.05 * size(a, 2));
a = a(mY+1:end-mY, mX+1:end-mX);
b = b(mY+1:end-mY, mX+1:end-mX);
a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
den = sqrt(sum(a.^2) * sum(b.^2));
if den == 0
    r = NaN;
else
    r = sum(a .* b) / den;
end

end

function tform = iMakeTform(A)
%IMAKETFORM Build a 2-D affine transform from a column-vector 3x3 matrix.

if exist('affinetform2d', 'class') == 8 || exist('affinetform2d', 'file') ~= 0
    tform = affinetform2d(A);
else
    tform = affine2d(A.');   % older releases use the row-vector convention
end

end

function p = iTformParams(tform)
%ITFORMPARAMS Return [tx ty rotationDeg scale] of a similarity transform.

if isprop(tform, 'A')
    A = tform.A;      % R2022b+: column-vector convention
else
    A = tform.T.';    % older releases: row-vector convention
end
p = [A(1,3), A(2,3), atan2d(A(2,1), A(1,1)), hypot(A(1,1), A(2,1))];

end

% =========================================================================
% Local helper: apply the estimated transforms to the .dat file
% =========================================================================
function [outPath, corrAfter] = iApplyTransformsToFile(entry, tforms, overwrite)
%IAPPLYTRANSFORMSTOFILE Warp every frame of the file with its transform.
%
% Also returns the correlation of each corrected frame with frame 1.

nt = entry.nt;
srcPath = entry.path;

if overwrite
    destPath = srcPath;
else
    [srcFolder, stem, ext] = fileparts(srcPath);
    destPath = fullfile(srcFolder, [stem '_MotionCorrected' ext]);
end

[destFolder, destStem, destExt] = fileparts(destPath);
tmpPath = fullfile(destFolder, [destStem '_writing' destExt]);
outputView = imref2d([entry.ny, entry.nx]);

% The corrected file is headered with the input's class, sizes, rate, and
% exposure, and the name of the file it is installed as.
try
    slabIn = spatialSlabIO('open', srcPath, 'Info', entry.info);
    cIn = onCleanup(@() spatialSlabIO('close', slabIn));
    slabOut = spatialSlabIO('create', tmpPath, datHeaderFromInfo(entry.info, destStem));
    cOut = onCleanup(@() spatialSlabIO('close', slabOut));

    corrAfter = ones(nt, 1);
    for t = 1:nt
        frame = iReadFrame(slabIn, entry, t);
        if t == 1
            fixed = iPrepareForRegistration(frame);
        end
        % NaNs are zeroed so a masked border behaves like the 0 fill of
        % pixels warped in from outside the frame.
        frame(isnan(frame)) = 0;
        if t > 1
            frame = imwarp(frame, tforms{t}, 'cubic', ...
                'OutputView', outputView, 'FillValues', 0);
            corrAfter(t) = iFrameCorrelation(fixed, iPrepareForRegistration(frame));
        end

        spatialSlabIO('write', slabOut, 1:entry.nx, frame, t);
    end
    spatialSlabIO('finalize', slabOut);
catch ME
    clear cIn cOut
    if isfile(tmpPath)
        delete(tmpPath);
    end
    rethrow(ME);
end

clear cIn cOut % close both handles before the file move below

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
% Local helpers: feedback and provenance
% =========================================================================
function iPrintSummary(entry, params, corrBefore, corrAfter, failedFrames)
%IPRINTSUMMARY Print the correction quality.

idx = 2:entry.nt;
nWorse = sum(corrAfter(idx) < corrBefore(idx) - 0.01);
fprintf(['%s: correlation with frame 1  mean %.3f -> %.3f | worst %.3f -> %.3f | ' ...
    'worse by >0.01 in %d/%d frames\n'], entry.name, ...
    mean(corrBefore(idx), 'omitnan'), mean(corrAfter(idx), 'omitnan'), ...
    min(corrBefore(idx)), min(corrAfter(idx)), nWorse, numel(idx));
fprintf('  max shift %.2f px | max rotation %.3f deg | scale %.4f-%.4f\n', ...
    max(hypot(params(:,1), params(:,2))), max(abs(params(:,3))), ...
    min(params(:,4)), max(params(:,4)));
if ~isempty(failedFrames)
    fprintf('  %d frame(s) could not be registered and were left uncorrected.\n', numel(failedFrames));
end

end

function iPlotQC(entry, transformType, params, corrBefore, corrAfter, showPlot, qcFile)
%IPLOTQC Plot correlation with frame 1 and the estimated motion.

isSimilarity = strcmp(transformType, 'similarity');
nCols = 2 + isSimilarity;
visible = 'off';
if showPlot
    visible = 'on';
end
fig = figure('Name', 'Motion correction QC', 'Visible', visible, ...
    'Position', [100 100 440*nCols 300]);
tl = tiledlayout(fig, 1, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s: motion correction (%s, cubic interpolation)', ...
    entry.name, transformType), 'Interpreter', 'none');

rate = NaN;
if isfield(entry.info, 'frameRateHz')
    rate = entry.info.frameRateHz;
end
if isnumeric(rate) && isscalar(rate) && isfinite(rate) && rate > 0
    x = (0:entry.nt-1) / rate;
    xLabel = 'Time (s)';
else
    x = 1:entry.nt;
    xLabel = 'Frame';
end

ax = nexttile(tl);
plot(ax, x, corrBefore, '-', 'Color', [0.6 0.6 0.6], 'DisplayName', 'Before'); hold(ax, 'on');
plot(ax, x, corrAfter, '-', 'Color', [0 0.447 0.741], 'DisplayName', 'After');
ylabel(ax, 'Correlation with frame 1'); xlabel(ax, xLabel); grid(ax, 'on');
legend(ax, 'Location', 'southwest');

ax = nexttile(tl);
plot(ax, x, params(:,2), '-', 'DisplayName', 'Row shift (Y)'); hold(ax, 'on');
plot(ax, x, params(:,1), '-', 'DisplayName', 'Column shift (X)');
ylabel(ax, 'Shift (px)'); xlabel(ax, xLabel); grid(ax, 'on');
legend(ax, 'Location', 'best');

if isSimilarity
    ax = nexttile(tl);
    plot(ax, x, params(:,3), '-');
    ylabel(ax, 'Rotation (deg)'); xlabel(ax, xLabel); grid(ax, 'on');
end

if ~isempty(qcFile)
    exportgraphics(fig, qcFile, 'Resolution', 120);
end
if ~showPlot
    close(fig);
end

end

function iSaveProvenance(SaveFolder, entry, transformType, usfac, params, ...
        corrBefore, corrAfter, failedFrames, outPath)
%ISAVEPROVENANCE Persist the transforms and run metadata for provenance/QC.

MotionCorrection = struct( ...
    'inputFile', entry.name, ...
    'transformType', transformType, ...
    'interpolation', 'cubic', ...
    'upsamplingFactor', usfac, ...
    'columns', {{'tx','ty','rotationDeg','scale'}}, ...
    'params', params, ...
    'corrBefore', corrBefore, ...
    'corrAfter', corrAfter, ...
    'failedFrames', failedFrames, ...
    'correctedFile', outPath, ...
    'appliedOn', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'appliedBy', mfilename);

save(fullfile(SaveFolder, [entry.stem '_MotionCorrection.mat']), 'MotionCorrection');

end

% =========================================================================
% Local pipeline info
% =========================================================================
function info = localPipelineInfo()
%LOCALPIPELINEINFO Return PipelineManager metadata for correctMotionArtifact.

info = PipelineManager.createPipelineInfo( ...
    'correctMotionArtifact', ...
    ['Estimate per-frame motion against the first frame of a .dat file ' ...
     'and correct it.']);

info = PipelineManager.addInput( ...
    info, ...
    'data', ...
    'ImageTimeSeries', ...
    'Y-X-T .dat filename to correct, estimated against its own first frame.', ...
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
    'Folder containing the .dat file; corrected files are written here.', ...
    'kind', 'input', ...
    'position', 2, ...
    'callType', 'positional', ...
    'isData', false);

info = PipelineManager.addInput( ...
    info, ...
    'TransformType', ...
    'parameter', ...
    'Motion model: translation (fast) or similarity (rotation, scale and translation).', ...
    'kind', 'parameter', ...
    'position', 3, ...
    'callType', 'namevalue', ...
    'default', 'translation', ...
    'allowed', {'translation','similarity'}, ...
    'dataType', 'char');

info = PipelineManager.addInput( ...
    info, ...
    'UpsamplingFactor', ...
    'parameter', ...
    'DFTREGISTRATION upsampling factor for translation (1 = whole-pixel, higher = subpixel).', ...
    'kind', 'parameter', ...
    'position', 4, ...
    'callType', 'namevalue', ...
    'default', 100, ...
    'allowed', [1 Inf], ...
    'dataType', 'double');

info = PipelineManager.addInput( ...
    info, ...
    'Overwrite', ...
    'parameter', ...
    ['If true, destructively rewrites the .dat file in place instead ' ...
     'of writing a "_MotionCorrected" copy.'], ...
    'kind', 'parameter', ...
    'position', 5, ...
    'callType', 'namevalue', ...
    'default', false, ...
    'allowed', [true false], ...
    'dataType', 'logical');

info = PipelineManager.addInput( ...
    info, ...
    'SaveShifts', ...
    'parameter', ...
    ['If true, saves the per-frame transforms and correlations to ' ...
     '<name>_MotionCorrection.mat and the QC figure to ' ...
     '<name>_MotionCorrectionQC.png.'], ...
    'kind', 'parameter', ...
    'position', 6, ...
    'callType', 'namevalue', ...
    'default', true, ...
    'allowed', [true false], ...
    'dataType', 'logical');

info = PipelineManager.addInput( ...
    info, ...
    'ShowPlot', ...
    'parameter', ...
    'If true, shows the correlation and motion QC figure.', ...
    'kind', 'parameter', ...
    'position', 7, ...
    'callType', 'namevalue', ...
    'default', true, ...
    'allowed', [true false], ...
    'dataType', 'logical');

info = PipelineManager.addOutput( ...
    info, ...
    'outData', ...
    'MotionShifts', ...
    'data', ...
    'Nt-by-4 [tx, ty, rotationDeg, scale] transform of each frame onto the first frame.', ...
    'MotionCorrectionShifts.mat', ...
    1, ...
    'isData', true, ...
    'saveFileName', 'MotionCorrectionShifts.mat');

info = PipelineManager.addOutput( ...
    info, ...
    'correctedDatFiles', ...
    {'ImageTimeSeries','ProcessedData'}, ...
    'file', ...
    ['The corrected .dat file, written as a new file ' ...
     '(or in place when Overwrite is true).'], ...
    '*.dat', ...
    2, ...
    'isData', false, ...
    'returnsValue', false);

info.notes = { ...
    ['Uses the third-party DFTREGISTRATION algorithm (Guizar-Sicairos, ' ...
     'Thurman & Fienup, Opt. Lett. 33, 156-158, 2008), bundled as a ' ...
     'private helper.']; ...
    'Overwrite=false (default) never modifies the source .dat file.'; ...
    'Each file is corrected independently against its own first frame.'};

end
