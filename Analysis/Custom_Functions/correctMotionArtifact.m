function outData = correctMotionArtifact(data, SaveFolder, varargin)
%CORRECTMOTIONARTIFACT Correct frame-wise motion artifacts in a YXT array or .dat file.
%
%   outData = correctMotionArtifact(data, SaveFolder)
%   outData = correctMotionArtifact(data, SaveFolder, Name, Value, ...)
%   info = correctMotionArtifact('pipelineInfo')
%
%   Estimates the transform that registers every frame of a Y-by-X-by-T
%   recording to its own first frame, and warps each frame with it using
%   cubic interpolation. Frames are high-pass filtered before registration
%   so that uneven illumination does not bias the estimate.
%
%   Each recording is corrected independently, against its own first frame.
%   The transforms are not shared between recordings: with strobed
%   illumination the channels are sampled at different instants, so motion
%   can differ between them. To correct several channels, run the function
%   once per channel.
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
%   data       - Either a numeric Y-X-T array with at least 2 frames, or the
%                filename of a continuous Y-X-T single-precision .dat file
%                (resolved relative to SaveFolder when not found as given).
%                Event-split (Y-X-T-E) data, UMT structs, and .umt files are
%                not supported.
%   SaveFolder - Existing folder. The corrected .dat file and the optional
%                run files are written here.
%
% Name-Value Options
%   TransformType    - 'translation' (default) or 'similarity'.
%   UpsamplingFactor - Positive integer upsampling factor of the
%                      'translation' estimate. 1 registers to the nearest
%                      whole pixel; higher values register to within
%                      1/factor of a pixel. Ignored for 'similarity'.
%                      Default: 100.
%   SaveShifts       - Logical scalar. When true (default), the shift matrix
%                      and the run metadata are saved to
%                      "<name>_MotionCorrection.mat" (variable
%                      MotionCorrection; its field "params" is the Nt-by-4
%                      matrix [tx, ty, rotationDeg, scale] of each frame
%                      onto the first frame, in pixels and degrees, with
%                      row 1 always [0 0 0 1] and rotation 0 and scale 1 for
%                      'translation'), and the QC figure to
%                      "<name>_MotionCorrectionQC.png", in SaveFolder. For
%                      array input the files are named "MotionCorrection.mat"
%                      and "MotionCorrectionQC.png". When false, neither
%                      file is written.
%   ShowPlot         - Logical scalar. When true (default), shows the QC
%                      figure.
%
% Output
%   outData - The corrected data, with the size of the input.
%             - Array input: a single-precision Y-X-T array.
%             - .dat input: the full path of "motionCorrected.dat" in
%               SaveFolder (same axes and class as the input). The source
%               file is never modified.
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
%   corrected = correctMotionArtifact('green.dat', saveFolder, ...
%       'TransformType', 'similarity');
%
% See also CREATEREGISTRATIONFORM, APPLYREGISTRATIONFORMONFOLDER.

default_Output = 'motionCorrected.dat';

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo(default_Output);
    return
end

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'TransformType', 'translation', ...
    @(x) (ischar(x) || (isstring(x) && isscalar(x))) && ...
    ismember(lower(char(string(x))), {'translation','similarity'}));
addParameter(p, 'UpsamplingFactor', 100, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == round(x));
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
saveShifts = p.Results.SaveShifts;
showPlot = p.Results.ShowPlot;

% Show warnings as a single line, without the call stack.
warnState = warning('off', 'backtrace');
restoreWarnings = onCleanup(@() warning(warnState));

entry = iResolveInput(data, SaveFolder);
if entry.nt < 2
    error('Umitoolbox:correctMotionArtifact:TooFewFrames', ...
        '%s must contain at least 2 frames.', entry.name);
end

fprintf('Estimating %s motion in %s...\n', transformType, entry.name);
[tforms, shifts, failedFrames, corrBefore] = iEstimateTransforms(entry, transformType, usfac);

fprintf('Correcting %s...\n', entry.name);
[outData, corrAfter] = iApplyTransforms(entry, tforms, SaveFolder, default_Output);

iPrintSummary(entry, shifts, corrBefore, corrAfter, failedFrames);

qcFile = '';
if saveShifts
    qcFile = fullfile(SaveFolder, [entry.sidecar 'QC.png']);
end
if saveShifts || showPlot
    iPlotQC(entry, transformType, shifts, corrBefore, corrAfter, showPlot, qcFile);
end

if saveShifts
    correctedFile = '';
    if strcmp(entry.kind, 'dat')
        correctedFile = outData;
    end
    iSaveProvenance(SaveFolder, entry, transformType, usfac, shifts, ...
        corrBefore, corrAfter, failedFrames, correctedFile);
end

end

% =========================================================================
% Local helpers: validate and resolve the input
% =========================================================================
function entry = iResolveInput(data, SaveFolder)
%IRESOLVEINPUT Resolve a Y-X-T array or a .dat filename into an entry struct.

if isnumeric(data) || islogical(data)
    validateattributes(data, {'numeric','logical'}, {'nonempty'}, mfilename, 'data');
    if ndims(data) ~= 3
        error('Umitoolbox:correctMotionArtifact:UnsupportedLayout', ...
            'Numeric input must be a Y x X x T array, got %d dimensions.', ndims(data));
    end
    entry = struct('kind', 'array', 'name', 'The input array', 'stem', '', ...
        'sidecar', 'MotionCorrection', 'path', '', ...
        'ny', size(data, 1), 'nx', size(data, 2), 'nt', size(data, 3), ...
        'info', struct(), 'data', single(data));
    return
end

if ~(ischar(data) || (isstring(data) && isscalar(data)))
    error('Umitoolbox:correctMotionArtifact:UnsupportedInputType', ...
        ['Input "data" must be a Y-X-T array or a .dat filename. ' ...
         'UMT structs are not supported.']);
end

inFile = char(string(data));
if ~isfile(inFile)
    altPath = fullfile(SaveFolder, inFile);
    if isfile(altPath)
        inFile = altPath;
    else
        error('Umitoolbox:correctMotionArtifact:InputFileNotFound', ...
            'Input file "%s" was not found.', char(string(data)));
    end
end

[~, inStem, inExt] = fileparts(inFile);
if ~strcmpi(inExt, '.dat')
    error('Umitoolbox:correctMotionArtifact:UnsupportedInputFile', ...
        'Unsupported input file extension "%s". Only .dat files are supported.', inExt);
end

entry = iResolveDatFile(inFile, [inStem inExt]);

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
entry = struct('kind', 'dat', 'name', fileName, 'stem', stem, ...
    'sidecar', [stem '_MotionCorrection'], 'path', datPath, ...
    'ny', datAxisSize(md, 'Y'), 'nx', datAxisSize(md, 'X'), ...
    'nt', datAxisSize(md, 'T'), 'info', md, 'data', []);

end

% =========================================================================
% Local helpers: estimate the transform of every frame onto frame 1
% =========================================================================
function [tforms, params, failedFrames, corrBefore] = iEstimateTransforms(entry, transformType, usfac)
%IESTIMATETRANSFORMS Estimate the transform of each frame onto frame 1.
%
% Also returns the correlation of each (uncorrected) frame with frame 1.

slabIn = [];
if strcmp(entry.kind, 'dat')
    slabIn = spatialSlabIO('open', entry.path, 'Info', entry.info);
    c = onCleanup(@() spatialSlabIO('close', slabIn));
end

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

if strcmp(entry.kind, 'array')
    frame = entry.data(:, :, t);
else
    frame = reshape(spatialSlabIO('read', slabIn, 1:entry.nx, t), entry.ny, entry.nx);
end

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
% Local helper: apply the estimated transforms
% =========================================================================
function [outData, corrAfter] = iApplyTransforms(entry, tforms, SaveFolder, defaultOutput)
%IAPPLYTRANSFORMS Warp every frame with its transform.
%
% Array input returns the corrected single array. A .dat input is written
% to DEFAULTOUTPUT in SaveFolder through a scratch file (so the declared
% output only appears once the run has completed, and the input may be the
% file it would overwrite) and its full path is returned. Also returns the
% correlation of each corrected frame with frame 1.

nt = entry.nt;
isDat = strcmp(entry.kind, 'dat');
outputView = imref2d([entry.ny, entry.nx]);
corrAfter = ones(nt, 1);

if isDat
    [~, outStem, outExt] = fileparts(defaultOutput);
    destPath = fullfile(SaveFolder, defaultOutput);
    tmpPath = fullfile(SaveFolder, [outStem '_writing' outExt]);
else
    outData = zeros(entry.ny, entry.nx, nt, 'single');
end

% The corrected file is headered with the input's class, sizes, rate, and
% exposure, and the name of the file it is installed as.
try
    if isDat
        slabIn = spatialSlabIO('open', entry.path, 'Info', entry.info);
        cIn = onCleanup(@() spatialSlabIO('close', slabIn));
        slabOut = spatialSlabIO('create', tmpPath, datHeaderFromInfo(entry.info, outStem));
        cOut = onCleanup(@() spatialSlabIO('close', slabOut));
    else
        slabIn = [];
    end

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

        if isDat
            spatialSlabIO('write', slabOut, 1:entry.nx, frame, t);
        else
            outData(:, :, t) = frame;
        end
    end

    if isDat
        spatialSlabIO('finalize', slabOut);
    end
catch ME
    clear cIn cOut
    if isDat && isfile(tmpPath)
        delete(tmpPath);
    end
    rethrow(ME);
end

if isDat
    clear cIn cOut % close both handles before the file move below

    [moveOk, moveMsg] = movefile(tmpPath, destPath, 'f');
    assert(moveOk, 'Umitoolbox:correctMotionArtifact:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', tmpPath, destPath, moveMsg);

    outData = destPath;
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

save(fullfile(SaveFolder, [entry.sidecar '.mat']), 'MotionCorrection');

end

% =========================================================================
% Local pipeline info
% =========================================================================
function info = localPipelineInfo(defaultOutput)
%LOCALPIPELINEINFO Return PipelineManager metadata for correctMotionArtifact.

info = PipelineManager.createPipelineInfo( ...
    'correctMotionArtifact', ...
    ['Estimate per-frame motion against the first frame of a Y-X-T array ' ...
     'or .dat file and correct it.']);

info = PipelineManager.addInput( ...
    info, ...
    'data', ...
    {'ImageTimeSeries','ProcessedData'}, ...
    'Y-X-T array or .dat filename to correct, estimated against its own first frame.', ...
    'kind', 'input', ...
    'position', 1, ...
    'callType', 'positional', ...
    'isData', true, ...
    'supportsFile', true, ...
    'dataMode', 'either');

info = PipelineManager.addInput( ...
    info, ...
    'SaveFolder', ...
    'SaveFolder', ...
    'Folder where the corrected .dat file and the optional run files are written.', ...
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
    'SaveShifts', ...
    'parameter', ...
    ['If true, saves the shift matrix (Nt-by-4 [tx, ty, rotationDeg, ' ...
     'scale]), correlations and run metadata to <name>_MotionCorrection.mat ' ...
     'and the QC figure to <name>_MotionCorrectionQC.png; nothing is ' ...
     'saved when false.'], ...
    'kind', 'parameter', ...
    'position', 5, ...
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
    'position', 6, ...
    'callType', 'namevalue', ...
    'default', true, ...
    'allowed', [true false], ...
    'dataType', 'logical');

info = PipelineManager.addOutput( ...
    info, ...
    'outData', ...
    {'ImageTimeSeries','ProcessedData'}, ...
    'data', ...
    'Motion-corrected Y-X-T data (array for array input, .dat for .dat input).', ...
    defaultOutput, ...
    1, ...
    'isData', true);

info.notes = { ...
    ['Uses the third-party DFTREGISTRATION algorithm (Guizar-Sicairos, ' ...
     'Thurman & Fienup, Opt. Lett. 33, 156-158, 2008), bundled as a ' ...
     'private helper.']; ...
    'The source .dat file is never modified.'; ...
    'Each recording is corrected independently against its own first frame.'};

end
