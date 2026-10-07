function outFile = correctMotionArtifact(SaveFolder, varargin)
%CORRECTMOTIONARTIFACT Correct frame-wise motion artifacts in imported movies.
%
%   outFile = correctMotionArtifact(SaveFolder)
%   outFile = correctMotionArtifact(SaveFolder, opts)
%
%   Independently for every imported Y-X-T .dat movie, estimates the
%   transform that registers each frame to the movie's own first frame, and
%   warps the frame with it using cubic interpolation. Frames are
%   high-pass filtered before registration so that uneven illumination does
%   not bias the estimate. Transforms are NOT shared between channels: with
%   strobed illumination the channels are sampled at different instants, so
%   motion can differ between them.
%
%   Two transform types are available (opts.TransformType):
%     'translation' - row/column shift estimated with the FFT-based
%                     DFTREGISTRATION algorithm (Guizar-Sicairos, Thurman &
%                     Fienup, Opt. Lett. 33, 156-158, 2008). Fast.
%     'similarity'  - rotation, isotropic scale and translation estimated by
%                     phase correlation (IMREGCORR). About 5x slower.
%
%   When finished, a summary of the correlation of every frame with the
%   first frame, before and after correction, is printed and plotted.
%
%   Imported channels are discovered from AcqInfos.mat (AcqInfoStream.
%   Illumination<N>.Color), as in TRIM_MOVIE. Every movie is processed in
%   RAM-safe mode, one frame at a time.
%
%   Inputs:
%       SaveFolder - Folder containing AcqInfos.mat, the imported .dat
%                    files and their same-name metadata .mat files.
%       opts       - Optional structure with fields:
%                    TransformType - 'translation' or 'similarity'.
%                                    Default: 'translation'
%                    Overwrite     - If false, each corrected movie is
%                                    written as "<name>_MotionCorrected.dat"
%                                    (with its own metadata .mat) and no
%                                    source file is modified. If true, every
%                                    channel .dat file is replaced in place.
%                                    Default: false
%                    SaveShifts    - If true, saves the per-frame transform
%                                    parameters, correlations and run info to
%                                    MotionCorrectionShifts.mat, and the QC
%                                    figure to MotionCorrectionQC.png, in
%                                    SaveFolder. Default: true
%                    ShowPlot      - If true, shows the QC figure.
%                                    Default: true
%
%   Output:
%       outFile - Cell array with the full paths of the corrected .dat
%                 files, one per imported channel.
%
%   Notes:
%       - Channels listed in AcqInfos.mat whose .dat or metadata .mat file is
%         missing are skipped with a warning; the remaining channels are
%         corrected. It is an error only if none of the files are found.
%       - All channel movies must be continuous Y-X-T, single precision,
%         with at least 2 frames. Every channel is validated before any file
%         is written.
%       - Cubic interpolation applies sub-pixel shifts without rounding, which
%         avoids the whole-pixel flips that nearest-neighbor interpolation
%         introduces between consecutive frames. It slightly smooths the
%         image, so corrected pixel values are interpolated, not original.
%       - Pixels warped in from outside the frame are filled with 0. NaN
%         pixels in the source are also set to 0 in the corrected movie.
%       - If the transform of a frame cannot be estimated, a warning is
%         issued and that frame is left uncorrected (identity transform).
%       - MotionCorrectionShifts.mat holds, per channel, an Nt-by-4 table
%         [tx, ty, rotationDeg, scale] in pixels/degrees (rotation 0 and
%         scale 1 for 'translation'). The first row is always [0 0 0 1].
%         Transforms map each frame onto the first frame.
%       - The reported correlation is the Pearson correlation between the
%         high-pass filtered first frame and each filtered frame, computed on
%         the image excluding a 5% border so that zero-filled edges do not
%         count against the correction. The summary counts frames whose
%         correlation dropped by more than 0.01; small drops on frames that
%         were already aligned come from rounding to whole pixels.
%
%   Example:
%       outFile = correctMotionArtifact(saveFolder, ...
%           struct('TransformType', 'similarity'));
%
%   See also TRIM_MOVIE, ALIGNFRAMES.

% Defaults
default_Output = {'fluo_475.dat', 'fluo_567.dat','fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'}; %#ok<NASGU> % Reference for PipelineManager. Actual outputs are stored in outFile.
default_opts = struct('TransformType', 'translation', 'Overwrite', false, 'SaveShifts', true, 'ShowPlot', true);
opts_values = struct('TransformType', {{'translation','similarity'}}, 'Overwrite', [true,false], 'SaveShifts', [true,false], 'ShowPlot', [true,false]);%#ok<NASGU>

%%% Arguments parsing and validation %%%
p = inputParser;
p.FunctionName = 'correctMotionArtifact';
addRequired(p, 'SaveFolder', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && isfolder(x));
addOptional(p, 'opts', default_opts, @(x) isstruct(x) && ~isempty(x));
parse(p, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
opts = p.Results.opts;
clear p

% Allow callers to provide partial opts structs.
optNames = fieldnames(default_opts);
for iOpt = 1:numel(optNames)
    if ~isfield(opts, optNames{iOpt}) || isempty(opts.(optNames{iOpt}))
        opts.(optNames{iOpt}) = default_opts.(optNames{iOpt});
    end
end

errID = 'umIToolbox:correctMotionArtifact:InvalidInput';
assert((ischar(opts.TransformType) || (isstring(opts.TransformType) && isscalar(opts.TransformType))) ...
    && ismember(lower(char(string(opts.TransformType))), {'translation','similarity'}), errID, ...
    'TransformType must be "translation" or "similarity".');
assert((islogical(opts.Overwrite) || isnumeric(opts.Overwrite)) && isscalar(opts.Overwrite), errID, ...
    'Overwrite must be a logical scalar.');
assert((islogical(opts.SaveShifts) || isnumeric(opts.SaveShifts)) && isscalar(opts.SaveShifts), errID, ...
    'SaveShifts must be a logical scalar.');
assert((islogical(opts.ShowPlot) || isnumeric(opts.ShowPlot)) && isscalar(opts.ShowPlot), errID, ...
    'ShowPlot must be a logical scalar.');
transformType = lower(char(string(opts.TransformType)));
overwrite = logical(opts.Overwrite);
saveShifts = logical(opts.SaveShifts);
showPlot = logical(opts.ShowPlot);

% Show warnings as a single line, without the call stack.
warnState = warning('off', 'backtrace');
restoreWarnings = onCleanup(@() warning(warnState));

%%% Discover and validate every channel before writing anything %%%
% A bad channel found partway through would leave others already corrected.
channelInfo = iDiscoverChannels(SaveFolder);
nChannels = numel(channelInfo);
for iChan = 1:nChannels
    assert(channelInfo(iChan).nt >= 2, 'umIToolbox:correctMotionArtifact:TooFewFrames', ...
        'File "%s" must contain at least 2 frames.', channelInfo(iChan).datFileName);
end

%%% Estimate and apply transforms, one channel at a time %%%
outFile = cell(1, nChannels);
channelParams = cell(1, nChannels);
failedFrames = cell(1, nChannels);
corrBefore = cell(1, nChannels);
corrAfter = cell(1, nChannels);
for iChan = 1:nChannels
    fprintf('Estimating %s motion in %s (%s)...\n', transformType, ...
        channelInfo(iChan).datFileName, channelInfo(iChan).colorName);
    [tforms, channelParams{iChan}, failedFrames{iChan}, corrBefore{iChan}] = ...
        iEstimateTransforms(channelInfo(iChan), transformType);
    fprintf('Correcting %s...\n', channelInfo(iChan).datFileName);
    [outFile{iChan}, corrAfter{iChan}] = ...
        iApplyTransformsToFile(channelInfo(iChan), tforms, overwrite);
    iPrintSummary(channelInfo(iChan), channelParams{iChan}, corrBefore{iChan}, ...
        corrAfter{iChan}, failedFrames{iChan});
end

if saveShifts || showPlot
    qcFile = '';
    if saveShifts
        qcFile = fullfile(SaveFolder, 'MotionCorrectionQC.png');
    end
    iPlotQC(channelInfo, transformType, channelParams, corrBefore, corrAfter, showPlot, qcFile);
end

if saveShifts
    MotionCorrection = struct( ...
        'files', {{channelInfo.datFileName}}, ...
        'transformType', transformType, ...
        'interpolation', 'cubic', ...
        'columns', {{'tx','ty','rotationDeg','scale'}}, ...
        'params', {channelParams}, ...
        'corrBefore', {corrBefore}, ...
        'corrAfter', {corrAfter}, ...
        'failedFrames', {failedFrames}, ...
        'correctedFiles', {outFile}, ...
        'appliedOn', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
        'appliedBy', mfilename);
    save(fullfile(SaveFolder, 'MotionCorrectionShifts.mat'), 'MotionCorrection');
end

disp('Done')
end

% =========================================================================
% Local helpers
% =========================================================================

function channelInfo = iDiscoverChannels(SaveFolder)
%IDISCOVERCHANNELS Resolve and validate every imported channel .dat file.

errID = 'umIToolbox:correctMotionArtifact:InvalidInput';

acqInfoFile = fullfile(SaveFolder, 'AcqInfos.mat');
assert(isfile(acqInfoFile), errID, ...
    'AcqInfos.mat was not found in SaveFolder: %s', SaveFolder);
acqInfo = load(acqInfoFile, 'AcqInfoStream', '-mat');
assert(isfield(acqInfo, 'AcqInfoStream') && isstruct(acqInfo.AcqInfoStream), errID, ...
    'AcqInfos.mat must contain the structure AcqInfoStream.');

illumFields = fieldnames(acqInfo.AcqInfoStream);
illumFields = illumFields(~cellfun(@isempty, regexp(illumFields, '^Illumination\d+$', 'once')));
assert(~isempty(illumFields), errID, ...
    'AcqInfoStream does not contain any Illumination<N> fields.');
illumIdx = cellfun(@(f) str2double(regexp(f, '\d+$', 'match', 'once')), illumFields);
[~, sortIdx] = sort(illumIdx);
illumFields = illumFields(sortIdx);

nChannels = numel(illumFields);
channelInfo = repmat(struct('colorName', '', 'datFileName', '', 'datFile', '', ...
    'matFile', '', 'metaData', struct(), 'ny', [], 'nx', [], 'nt', []), nChannels, 1);

for iChan = 1:nChannels
    illumInfo = acqInfo.AcqInfoStream.(illumFields{iChan});
    assert(isstruct(illumInfo) && isfield(illumInfo, 'Color'), errID, ...
        'AcqInfoStream.%s must contain the field Color.', illumFields{iChan});
    colorName = char(string(illumInfo.Color));
    datFileName = iColorNameToDatFileName(colorName);

    channelInfo(iChan).colorName = colorName;
    channelInfo(iChan).datFileName = datFileName;
    channelInfo(iChan).datFile = fullfile(SaveFolder, datFileName);
    channelInfo(iChan).matFile = fullfile(SaveFolder, strrep(datFileName, '.dat', '.mat'));
end

datFileNames = {channelInfo.datFileName};
for iChan = 1:nChannels
    assert(nnz(strcmp(datFileNames{iChan}, datFileNames)) == 1, errID, ...
        'Multiple illumination channels map to the same .dat file: %s.', datFileNames{iChan});
end

% Channels listed in AcqInfos.mat whose files were removed are skipped with
% a warning; it is an error only if no channel has its files.
isPresent = true(nChannels, 1);
for iChan = 1:nChannels
    missing = {};
    if ~isfile(channelInfo(iChan).datFile)
        missing{end+1} = channelInfo(iChan).datFileName; %#ok<AGROW>
    end
    if ~isfile(channelInfo(iChan).matFile)
        [~, matName, matExt] = fileparts(channelInfo(iChan).matFile);
        missing{end+1} = [matName matExt]; %#ok<AGROW>
    end
    if ~isempty(missing)
        isPresent(iChan) = false;
        warning('umIToolbox:correctMotionArtifact:MissingChannelFile', ...
            'Skipped %s: %s not found.', channelInfo(iChan).colorName, strjoin(missing, ', '));
    end
end
assert(any(isPresent), 'umIToolbox:correctMotionArtifact:NoFiles', ...
    'None of the imported channel files listed in AcqInfos.mat were found in "%s".', SaveFolder);
channelInfo = channelInfo(isPresent);
nChannels = numel(channelInfo);

for iChan = 1:nChannels
    md = load(channelInfo(iChan).matFile, '-mat');
    requiredFields = {'dim_names', 'datSize', 'datLength', 'Datatype'};
    missingFields = requiredFields(~isfield(md, requiredFields));
    assert(isempty(missingFields), errID, ...
        'Metadata file %s is missing required field(s): %s.', ...
        channelInfo(iChan).matFile, strjoin(missingFields, ', '));

    assert(isequal(cellstr(md.dim_names(:)).', {'Y','X','T'}), ...
        'umIToolbox:correctMotionArtifact:UnsupportedLayout', ...
        'File "%s" has dimensions {%s}. Only continuous Y-X-T movies are supported.', ...
        channelInfo(iChan).datFileName, strjoin(cellstr(md.dim_names(:)).', ','));
    assert(strcmp(md.Datatype, 'single'), ...
        'umIToolbox:correctMotionArtifact:UnsupportedDataClass', ...
        'Only single-precision movies are supported; "%s" stores %s.', ...
        channelInfo(iChan).datFileName, md.Datatype);

    channelInfo(iChan).metaData = md;
    channelInfo(iChan).ny = md.datSize(1);
    channelInfo(iChan).nx = md.datSize(2);
    channelInfo(iChan).nt = md.datLength(1);

    % The file size must match the metadata, or frame offsets are wrong.
    expectedBytes = 4 * double(channelInfo(iChan).ny) * channelInfo(iChan).nx * channelInfo(iChan).nt;
    fileInfo = dir(channelInfo(iChan).datFile);
    assert(fileInfo.bytes == expectedBytes, errID, ...
        '"%s" is %d bytes, but its metadata describes %d bytes.', ...
        channelInfo(iChan).datFileName, fileInfo.bytes, expectedBytes);
end

end

function datFileName = iColorNameToDatFileName(colorName)
%ICOLORNAMETODATFILENAME Convert an illumination color name to a .dat name.

colorName = strtrim(char(string(colorName)));
switch lower(colorName)
    case 'red'
        datFileName = 'red.dat';
    case {'amber', 'yellow'}
        datFileName = 'yellow.dat';
    case 'green'
        datFileName = 'green.dat';
    case 'fluo'
        datFileName = 'fluo.dat';
    case 'speckle'
        datFileName = 'speckle.dat';
    otherwise
        tokens = regexp(colorName, '^Fluo\s*#\d+\s+(\d+)\s*nm$', 'tokens', 'once');
        if isempty(tokens)
            error('umIToolbox:correctMotionArtifact:InvalidChannelName', ...
                'Unsupported illumination Color name: "%s".', colorName);
        end
        datFileName = sprintf('fluo_%s.dat', tokens{1});
end

end

function [tforms, params, failedFrames, corrBefore] = iEstimateTransforms(entry, transformType)
%IESTIMATETRANSFORMS Estimate the transform of each frame onto frame 1.
%
% Also returns the correlation of each (uncorrected) frame with frame 1.

fid = fopen(entry.datFile, 'rb');
assert(fid > 0, 'umIToolbox:correctMotionArtifact:FileOpenFailed', ...
    'Failed to open input file: %s', entry.datFile);
c = onCleanup(@() safeFclose(fid));

fixed = iPrepareForRegistration(iReadFrame(fid, entry, 1));
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
    moving = iPrepareForRegistration(iReadFrame(fid, entry, t));
    corrBefore(t) = iFrameCorrelation(fixed, moving);
    try
        switch transformType
            case 'translation'
                reg = dftregistration(fixedFFT, fft2(moving), 100);
                tx = reg(4);   % column shift
                ty = reg(3);   % row shift
                tforms{t} = iMakeTform([1 0 tx; 0 1 ty; 0 0 1]);
                params(t, :) = [tx ty 0 1];
            case 'similarity'
                tforms{t} = imregcorr(moving, fixed, 'similarity');
                params(t, :) = iTformParams(tforms{t});
        end
    catch ME
        warning('umIToolbox:correctMotionArtifact:RegistrationFailed', ...
            'Frame %d of %s was left uncorrected: %s', t, entry.datFileName, ME.message);
        tforms{t} = iMakeTform(eye(3));
        params(t, :) = [0 0 0 1];
        failedFrames(end+1) = t; %#ok<AGROW>
    end
    if mod(t, 500) == 0
        fprintf('  %s: %d/%d frames\n', entry.datFileName, t, entry.nt);
    end
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

function p = iTformParams(tform)
%ITFORMPARAMS Return [tx ty rotationDeg scale] of a similarity transform.

if isprop(tform, 'A')
    A = tform.A;      % R2022b+: column-vector convention
else
    A = tform.T.';    % older releases: row-vector convention
end
p = [A(1,3), A(2,3), atan2d(A(2,1), A(1,1)), hypot(A(1,1), A(2,1))];

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

function frame = iReadFrame(fid, entry, t)
%IREADFRAME Read frame t (1-based) as a Y-by-X single array.

nPix = double(entry.ny) * entry.nx;
status = fseek(fid, (t-1) * nPix * 4, 'bof');
assert(status == 0, 'umIToolbox:correctMotionArtifact:ReadFailed', ...
    'Failed to seek to frame %d in %s.', t, entry.datFile);
frame = fread(fid, nPix, '*single');
assert(numel(frame) == nPix, 'umIToolbox:correctMotionArtifact:ReadFailed', ...
    'Failed to read frame %d from %s.', t, entry.datFile);
frame = reshape(frame, entry.ny, entry.nx);

end

function [outPath, corrAfter] = iApplyTransformsToFile(entry, tforms, overwrite)
%IAPPLYTRANSFORMSTOFILE Warp every frame of one movie with its transform.
%
% Also returns the correlation of each corrected frame with frame 1.

if overwrite
    destPath = entry.datFile;
else
    [folder, stem, ext] = fileparts(entry.datFile);
    destPath = fullfile(folder, [stem '_MotionCorrected' ext]);
end
[destFolder, destStem, destExt] = fileparts(destPath);
tmpPath = fullfile(destFolder, [destStem '_writing' destExt]);
outputView = imref2d([entry.ny, entry.nx]);

try
    fidIn = fopen(entry.datFile, 'rb');
    assert(fidIn > 0, 'umIToolbox:correctMotionArtifact:FileOpenFailed', ...
        'Failed to open input file: %s', entry.datFile);
    cIn = onCleanup(@() safeFclose(fidIn));
    fidOut = fopen(tmpPath, 'wb');
    assert(fidOut > 0, 'umIToolbox:correctMotionArtifact:FileOpenFailed', ...
        'Failed to create temporary output file: %s', tmpPath);
    cOut = onCleanup(@() safeFclose(fidOut));

    corrAfter = ones(entry.nt, 1);
    for t = 1:entry.nt
        frame = iReadFrame(fidIn, entry, t);
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
        nWritten = fwrite(fidOut, frame, 'single');
        assert(nWritten == numel(frame), 'umIToolbox:correctMotionArtifact:WriteFailed', ...
            'Failed to write frame %d to %s.', t, tmpPath);
    end
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
    [ok, msg] = movefile(tmpPath, destPath, 'f');
    assert(ok, 'umIToolbox:correctMotionArtifact:OutputMoveFailed', ...
        'Failed to move "%s" onto "%s": %s', tmpPath, destPath, msg);
    metaData = entry.metaData;
    metaData.datFile = destPath;
    save(strrep(destPath, '.dat', '.mat'), '-struct', 'metaData');
end

outPath = destPath;

end

function iReplaceFileSafely(tmpPath, destPath)
%IREPLACEFILESAFELY Install tmpPath as destPath without a data-loss window.

backupPath = [destPath '.bak'];
backupCreated = false;

if isfile(destPath)
    [ok, msg] = movefile(destPath, backupPath, 'f');
    assert(ok, 'umIToolbox:correctMotionArtifact:BackupFailed', ...
        'Could not back up "%s" before replacing it: %s', destPath, msg);
    backupCreated = true;
end

[ok, msg] = movefile(tmpPath, destPath, 'f');
if ~ok
    if backupCreated && isfile(backupPath)
        movefile(backupPath, destPath, 'f');
    end
    error('umIToolbox:correctMotionArtifact:ReplaceFailed', ...
        'Could not install the corrected data as "%s": %s', destPath, msg);
end

if backupCreated && isfile(backupPath)
    delete(backupPath);
end

end

% =========================================================================
% Feedback to the user
% =========================================================================

function iPrintSummary(entry, params, corrBefore, corrAfter, failedFrames)
%IPRINTSUMMARY Print the correction quality of one channel.

idx = 2:entry.nt;
nWorse = sum(corrAfter(idx) < corrBefore(idx) - 0.01);
fprintf(['%s (%s): correlation with frame 1  mean %.3f -> %.3f | worst %.3f -> %.3f | ' ...
    'worse by >0.01 in %d/%d frames\n'], entry.datFileName, entry.colorName, ...
    mean(corrBefore(idx), 'omitnan'), mean(corrAfter(idx), 'omitnan'), ...
    min(corrBefore(idx)), min(corrAfter(idx)), nWorse, numel(idx));
fprintf('  max shift %.2f px | max rotation %.3f deg | scale %.4f-%.4f\n', ...
    max(hypot(params(:,1), params(:,2))), max(abs(params(:,3))), ...
    min(params(:,4)), max(params(:,4)));
if ~isempty(failedFrames)
    fprintf('  %d frame(s) could not be registered and were left uncorrected.\n', numel(failedFrames));
end

end

function iPlotQC(channelInfo, transformType, channelParams, corrBefore, corrAfter, showPlot, qcFile)
%IPLOTQC Plot correlation with frame 1 and the estimated motion per channel.

nChannels = numel(channelInfo);
isSimilarity = strcmp(transformType, 'similarity');
nCols = 2 + isSimilarity;
visible = 'off';
if showPlot
    visible = 'on';
end
fig = figure('Name', 'Motion correction QC', 'Visible', visible, ...
    'Position', [100 100 440*nCols 260*nChannels]);
tl = tiledlayout(fig, nChannels, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Motion correction (%s, cubic interpolation)', transformType));

for iChan = 1:nChannels
    entry = channelInfo(iChan);
    if isfield(entry.metaData, 'Freq') && isnumeric(entry.metaData.Freq) && ...
            isscalar(entry.metaData.Freq) && entry.metaData.Freq > 0
        x = (0:entry.nt-1) / entry.metaData.Freq;
        xLabel = 'Time (s)';
    else
        x = 1:entry.nt;
        xLabel = 'Frame';
    end
    params = channelParams{iChan};

    ax = nexttile(tl);
    plot(ax, x, corrBefore{iChan}, '-', 'Color', [0.6 0.6 0.6], 'DisplayName', 'Before'); hold(ax, 'on');
    plot(ax, x, corrAfter{iChan}, '-', 'Color', [0 0.447 0.741], 'DisplayName', 'After');
    ylabel(ax, 'Correlation with frame 1'); xlabel(ax, xLabel); grid(ax, 'on');
    title(ax, sprintf('%s (%s)', entry.datFileName, entry.colorName), 'Interpreter', 'none');
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
end

if ~isempty(qcFile)
    exportgraphics(fig, qcFile, 'Resolution', 120);
end
if ~showPlot
    close(fig);
end

end
