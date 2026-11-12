function outData = correctMotionArtifact(data, SaveFolder, varargin)
%CORRECTMOTIONARTIFACT Correct frame-wise motion artifacts in a YXT array or .dat file(s).
%
%   outData = correctMotionArtifact(data, SaveFolder)
%   outData = correctMotionArtifact(data, SaveFolder, Name, Value, ...)
%   outFiles = correctMotionArtifact(SaveFolder)
%   outFiles = correctMotionArtifact(SaveFolder, Name, Value, ...)
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
%   can differ between them. To correct several channels, either call the
%   function once per channel or use the batch form, which takes a folder
%   instead of data and corrects every .dat file listed in AcqInfos.mat
%   (see Batch mode).
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
%   data       - One of:
%                - a numeric Y-X-T array with at least 2 frames;
%                - the filename of a continuous Y-X-T single-precision .dat
%                  file (resolved relative to SaveFolder when not found as
%                  given);
%                Event-split (Y-X-T-E) data, UMT structs, and .umt files are
%                not supported. In the batch form this argument is left out
%                (see Batch mode).
%   SaveFolder - Existing folder. The corrected .dat file and the optional
%                run files are written here. In batch mode it is also the
%                folder that holds AcqInfos.mat and the channel files.
%
% Name-Value Options
%   TransformType    - 'translation' (default) or 'similarity'.
%   UpsamplingFactor - Positive integer upsampling factor of the
%                      'translation' estimate. 1 registers to the nearest
%                      whole pixel; higher values register to within
%                      1/factor of a pixel. Ignored for 'similarity'.
%                      Default: 100.
%   Overwrite        - Logical scalar. When false (default), each corrected
%                      .dat file is written to SaveFolder as
%                      "<name>_MotionCorrected.dat" and no source file is
%                      modified. When true, the source .dat file is replaced
%                      by its corrected version in place (nothing else is
%                      renamed or removed). Ignored for array input.
%   AllowRecorrection - Logical scalar. By default (false) a .dat file that
%                      was already motion corrected is refused, because
%                      correcting it again compounds the smoothing of the
%                      cubic interpolation without recovering any motion
%                      (see Already-corrected guard). When true, such a file
%                      is processed anyway after a warning. Default: false.
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
%             - .dat input: the full path of the corrected .dat file (same
%               axes and class as the input): "<name>_MotionCorrected.dat" in
%               SaveFolder, or the source file itself when Overwrite is true.
%             - Batch mode: a cell array with one such full path per
%               corrected channel.
%
% Batch mode
%   outFiles = correctMotionArtifact(SaveFolder, ...) is called with the
%   folder as its first argument (a .dat file is never a folder, so the form
%   is unambiguous). It corrects, one after the other, every .dat file
%   listed in the ImportedChannels of
%   SaveFolder's AcqInfos.mat, each against its own first frame and with its
%   own transforms, using the same options for all of them. Before any file
%   is written, every channel is checked (continuous Y-X-T, single
%   precision, at least 2 frames), so a bad channel stops the run with
%   nothing modified. A listed channel whose .dat file is missing is skipped
%   with a warning, and so is a channel that was already motion corrected
%   (unless AllowRecorrection is true); the call is an error only when no
%   channel is left to correct.
%
%   ImportedChannels is authoritative: only the files it lists are
%   corrected. A legacy folder whose AcqInfos.mat has no ImportedChannels
%   has its channels taken from the declared Illumination<N> entries, named
%   as the importer names them: red.dat, green.dat, yellow.dat (amber),
%   speckle.dat, fluo_<nm>.dat (a Fluo color with a 3-digit wavelength) and
%   fluo.dat. When nothing is declared, the red, green, yellow, speckle and
%   fluo files found in the folder are used. A declared channel whose file
%   is missing is skipped with a warning.
%
%   Batch mode is the form PipelineManager runs. The function is a pipeline
%   setup step: it has no data input or output, runs before the processing
%   steps, and with Overwrite=true it keeps the file names that
%   ImportedChannels lists, so every later step (HemoCompute, GSR, ...)
%   reads the motion-corrected files. With Overwrite=false it writes the
%   "_MotionCorrected" copies, which later steps do not use. The array and
%   single-file forms are for direct calls only.
%
% Already-corrected guard
%   A .dat file counts as already motion corrected when its name ends in
%   "_MotionCorrected", or when it is a file this function wrote earlier and
%   has not changed since, including a source file replaced in place with
%   Overwrite=true, which has no suffix to tell it apart. The function
%   records every .dat file it writes (name, size and modification time) in
%   "MotionCorrectionLog.mat", next to that file; the record is kept even
%   when SaveShifts is false. Replacing or regenerating the file, for
%   example by importing it again, makes it a new file that is no longer
%   recognized as corrected. Copies renamed outside this function, and files
%   corrected by older versions that kept no record, cannot be recognized.
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
%   - The corrected file is written to a scratch file and moved into place
%     only once complete, so a failed run never leaves a partial file, and
%     an in-place replacement is never half done.
%   - The reported correlation is the Pearson correlation between the
%     high-pass filtered first frame and each filtered frame, computed on
%     the image excluding a 5% border so that zero-filled edges do not count
%     against the correction. The summary counts frames whose correlation
%     dropped by more than 0.01.
%   - Overwrite=true is irreversible: the uncorrected recording is gone
%     (the guard only prevents correcting it a second time). Keep the raw
%     data to import again.
%
% Example
%   corrected = correctMotionArtifact('green.dat', saveFolder, ...
%       'TransformType', 'similarity');
%   allFiles = correctMotionArtifact(saveFolder);
%   correctMotionArtifact(saveFolder, 'Overwrite', true);
%
% See also CREATEREGISTRATIONFORM, APPLYREGISTRATIONFORMONFOLDER.

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

% Batch form: the folder comes first and there is no data argument.
isBatch = false;
if nargin >= 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && isfolder(char(string(data)))
    isBatch = true;
    if nargin >= 2
        varargin = [{SaveFolder}, varargin];
    end
    SaveFolder = data;
    data = [];
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
addParameter(p, 'Overwrite', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'AllowRecorrection', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'SaveShifts', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ShowPlot', true, @(x) islogical(x) && isscalar(x));
parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
if ~isfolder(SaveFolder)
    error('Umitoolbox:correctMotionArtifact:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

opts = struct( ...
    'transformType', lower(char(string(p.Results.TransformType))), ...
    'usfac', double(p.Results.UpsamplingFactor), ...
    'overwrite', p.Results.Overwrite, ...
    'allowRecorrection', p.Results.AllowRecorrection, ...
    'saveShifts', p.Results.SaveShifts, ...
    'showPlot', p.Results.ShowPlot);

% Show warnings as a single line, without the call stack.
warnState = warning('off', 'backtrace');
restoreWarnings = onCleanup(@() warning(warnState));

if isBatch
    outData = iCorrectFolder(SaveFolder, opts);
    return
end

entry = iResolveInput(data, SaveFolder);
iValidateFrameCount(entry);
reason = iPriorCorrection(entry);
if ~isempty(reason)
    if ~opts.allowRecorrection
        error('Umitoolbox:correctMotionArtifact:AlreadyCorrected', ...
            ['%s was already motion corrected (%s). Correcting it again ' ...
             'only blurs it further. Use the original recording, or set ' ...
             '''AllowRecorrection'' to true to proceed anyway.'], entry.name, reason);
    end
    warning('Umitoolbox:correctMotionArtifact:AlreadyCorrected', ...
        '%s was already motion corrected (%s); correcting it again.', entry.name, reason);
end

outData = iCorrectEntry(entry, SaveFolder, opts);

end

% =========================================================================
% Local helpers: run one recording, run a folder
% =========================================================================
function outData = iCorrectEntry(entry, SaveFolder, opts)
%ICORRECTENTRY Estimate and apply the correction of one recording, then
% save its QC and provenance files.

fprintf('Estimating %s motion in %s...\n', opts.transformType, entry.name);
[tforms, shifts, failedFrames, corrBefore] = ...
    iEstimateTransforms(entry, opts.transformType, opts.usfac);

destPath = '';
if strcmp(entry.kind, 'dat')
    if opts.overwrite
        destPath = entry.path;
    else
        destPath = fullfile(SaveFolder, [entry.stem '_MotionCorrected.dat']);
    end
end

fprintf('Correcting %s...\n', entry.name);
[outData, corrAfter] = iApplyTransforms(entry, tforms, destPath);
if strcmp(entry.kind, 'dat')
    iRecordCorrection(destPath, entry, opts.transformType);
end

iPrintSummary(entry, shifts, corrBefore, corrAfter, failedFrames);

qcFile = '';
if opts.saveShifts
    qcFile = fullfile(SaveFolder, [entry.sidecar 'QC.png']);
end
if opts.saveShifts || opts.showPlot
    iPlotQC(entry, opts.transformType, shifts, corrBefore, corrAfter, opts.showPlot, qcFile);
end

if opts.saveShifts
    iSaveProvenance(SaveFolder, entry, opts.transformType, opts.usfac, shifts, ...
        corrBefore, corrAfter, failedFrames, destPath);
end

end

function outFiles = iCorrectFolder(SaveFolder, opts)
%ICORRECTFOLDER Correct every .dat file listed in the folder's AcqInfos.mat.
%
% Every channel is resolved and checked before any file is written.

datNames = iListImportedDatFiles(SaveFolder);

entries = {};
nSkippedCorrected = 0;
for k = 1:numel(datNames)
    datPath = fullfile(SaveFolder, datNames{k});
    if ~isfile(datPath)
        warning('Umitoolbox:correctMotionArtifact:MissingChannelFile', ...
            'Skipped %s: the file was not found in "%s".', datNames{k}, SaveFolder);
        continue
    end
    entry = iResolveDatFile(datPath, datNames{k});
    iValidateFrameCount(entry);
    reason = iPriorCorrection(entry);
    if ~isempty(reason)
        if ~opts.allowRecorrection
            warning('Umitoolbox:correctMotionArtifact:AlreadyCorrected', ...
                'Skipped %s: it was already motion corrected (%s).', entry.name, reason);
            nSkippedCorrected = nSkippedCorrected + 1;
            continue
        end
        warning('Umitoolbox:correctMotionArtifact:AlreadyCorrected', ...
            '%s was already motion corrected (%s); correcting it again.', entry.name, reason);
    end
    entries{end+1} = entry; %#ok<AGROW>
end

if isempty(entries)
    if nSkippedCorrected > 0
        error('Umitoolbox:correctMotionArtifact:AlreadyCorrected', ...
            ['Every channel listed in AcqInfos.mat was already motion corrected. ' ...
             'Set ''AllowRecorrection'' to true to correct them again.']);
    end
    error('Umitoolbox:correctMotionArtifact:NoFiles', ...
        'None of the channel files listed in AcqInfos.mat were found in "%s".', SaveFolder);
end

outFiles = cell(1, numel(entries));
for k = 1:numel(entries)
    outFiles{k} = iCorrectEntry(entries{k}, SaveFolder, opts);
end

end

function datNames = iListImportedDatFiles(SaveFolder)
%ILISTIMPORTEDDATFILES Names of the imported .dat files of a folder.
%
% AcqInfoStream.ImportedChannels is authoritative when it exists: it lists
% exactly the files the importer wrote, so a derived file that merely sits
% in the folder is never touched. A folder that predates that field has the
% channels inferred by iInferLegacyDatFiles.

errID = 'Umitoolbox:correctMotionArtifact:InvalidInput';
acqFile = fullfile(SaveFolder, 'AcqInfos.mat');
if ~isfile(acqFile)
    error(errID, 'AcqInfos.mat was not found in SaveFolder: %s', SaveFolder);
end
S = load(acqFile, 'AcqInfoStream');
if ~isfield(S, 'AcqInfoStream') || ~isstruct(S.AcqInfoStream) || ~isscalar(S.AcqInfoStream)
    error(errID, 'AcqInfos.mat must contain the scalar structure AcqInfoStream.');
end

acqInfo = S.AcqInfoStream;
if isfield(acqInfo, 'ImportedChannels') && ~isempty(acqInfo.ImportedChannels) ...
        && isfield(acqInfo.ImportedChannels, 'DatFile')
    datNames = cell(1, numel(acqInfo.ImportedChannels));
    for k = 1:numel(datNames)
        datNames{k} = iDatFileName(acqInfo.ImportedChannels(k).DatFile);
    end
else
    datNames = iInferLegacyDatFiles(acqInfo, SaveFolder);
end
datNames = unique(datNames(~cellfun(@isempty, datNames)), 'stable');

if isempty(datNames)
    error(errID, 'AcqInfos.mat in "%s" does not list any imported .dat file.', SaveFolder);
end

end

function datNames = iInferLegacyDatFiles(acqInfo, SaveFolder)
%IINFERLEGACYDATFILES Channel files of a folder without ImportedChannels.
%
% Collects, in this order, and merges:
%   1. what resolveImportedChannelFallback infers (red, green, yellow). It
%      stops with an error when its evidence is incomplete, which here only
%      means less evidence, so its error is ignored;
%   2. every declared Illumination<N>, named the way the importer names its
%      files (red, green, yellow for amber, speckle, fluo_<nm>, fluo), which
%      covers the channels the fallback does not know;
%   3. only when nothing is declared, the canonical files found in the folder.
% A name whose file does not exist is reported by the caller, with a warning.

datNames = {};

try
    resolved = resolveImportedChannelFallback(acqInfo, SaveFolder);
    if isfield(resolved, 'ImportedChannels') && ~isempty(resolved.ImportedChannels)
        for k = 1:numel(resolved.ImportedChannels)
            datNames{end+1} = iDatFileName(resolved.ImportedChannels(k).DatFile); %#ok<AGROW>
        end
    end
catch
    % Incomplete or ambiguous evidence for the channels the fallback knows.
end

fields = fieldnames(acqInfo);
tokens = regexp(fields, '^Illumination(\d+)$', 'tokens', 'once');
isIllum = ~cellfun(@isempty, tokens);
fields = fields(isIllum);
[~, order] = sort(cellfun(@(t) str2double(t{1}), tokens(isIllum)));
for k = order(:).'
    illum = acqInfo.(fields{k});
    if ~isstruct(illum) || ~isscalar(illum)
        continue
    end
    if isfield(illum, 'DatFile') && ~isempty(illum.DatFile)
        datNames{end+1} = iDatFileName(illum.DatFile); %#ok<AGROW>
    elseif isfield(illum, 'Color') && ~isempty(illum.Color)
        datNames{end+1} = [iChannelTag(illum.Color) '.dat']; %#ok<AGROW>
    end
end

if isempty(datNames)
    for tag = {'red', 'green', 'yellow', 'speckle', 'fluo'}
        if isfile(fullfile(SaveFolder, [tag{1} '.dat']))
            datNames{end+1} = [tag{1} '.dat']; %#ok<AGROW>
        end
    end
    fluoFiles = dir(fullfile(SaveFolder, 'fluo_*.dat'));
    datNames = [datNames, {fluoFiles.name}];
end

end

function tag = iChannelTag(colorName)
%ICHANNELTAG Channel file tag of an illumination colour, as the importer
% (ImagesClassification) names it.

colorName = char(string(colorName));
if contains(colorName, 'red', 'IgnoreCase', true)
    tag = 'red';
elseif contains(colorName, 'green', 'IgnoreCase', true)
    tag = 'green';
elseif contains(colorName, {'amber', 'yellow'}, 'IgnoreCase', true)
    tag = 'yellow';
elseif contains(colorName, 'speckle', 'IgnoreCase', true)
    tag = 'speckle';
else
    waveTag = regexp(colorName, '[0-9]{3}', 'match');
    if ~isempty(waveTag)
        tag = ['fluo_' waveTag{:}];
    else
        tag = 'fluo';
    end
end

end

function name = iDatFileName(datFile)
%IDATFILENAME File name with extension (.dat when none) of a listed file.

[~, base, ext] = fileparts(char(string(datFile)));
if isempty(ext)
    ext = '.dat';
end
name = [base ext];

end


function iValidateFrameCount(entry)
%IVALIDATEFRAMECOUNT Registration needs a reference frame and at least one more.

if entry.nt < 2
    error('Umitoolbox:correctMotionArtifact:TooFewFrames', ...
        '%s must contain at least 2 frames.', entry.name);
end

end

% =========================================================================
% Local helpers: the already-corrected guard
% =========================================================================
function reason = iPriorCorrection(entry)
%IPRIORCORRECTION Explain why a .dat file counts as already corrected.
%
% Returns '' when it does not. A file counts as corrected when its name
% carries the "_MotionCorrected" suffix, or when MotionCorrectionLog.mat
% (next to the file) lists a file of this name whose size and modification
% time are unchanged since this function wrote it.

reason = '';
if ~strcmp(entry.kind, 'dat')
    return
end

if endsWith(entry.stem, '_MotionCorrected', 'IgnoreCase', true)
    reason = 'its name ends in "_MotionCorrected"';
    return
end

[folder, name, ext] = fileparts(entry.path);
logFile = fullfile(folder, 'MotionCorrectionLog.mat');
if ~isfile(logFile)
    return
end
try
    S = load(logFile, 'MotionCorrectionLog');
    log = S.MotionCorrectionLog;
catch
    warning('Umitoolbox:correctMotionArtifact:UnreadableLog', ...
        'Could not read "%s"; the check for earlier corrections was skipped.', logFile);
    return
end

d = dir(entry.path);
for k = 1:numel(log)
    if strcmpi(log(k).file, [name ext]) && log(k).bytes == d.bytes ...
            && abs(log(k).modified - iModifiedMs(entry.path)) < 1
        reason = sprintf('corrected on %s by this function, and unchanged since', ...
            log(k).appliedOn);
        return
    end
end

end

function ms = iModifiedMs(filePath)
%IMODIFIEDMS Modification time of a file in milliseconds.
%
% DIR reports the modification time to the whole second only, which cannot
% tell a file from one regenerated a moment later.

ms = double(java.io.File(filePath).lastModified());

end

function iRecordCorrection(destPath, entry, transformType)
%IRECORDCORRECTION Record a .dat file this function just wrote.
%
% The record lives in MotionCorrectionLog.mat in the file's folder and is
% what lets iPriorCorrection recognize a file corrected in place.

[folder, name, ext] = fileparts(destPath);
logFile = fullfile(folder, 'MotionCorrectionLog.mat');
d = dir(destPath);
rec = struct('file', [name ext], 'bytes', d.bytes, 'modified', iModifiedMs(destPath), ...
    'inputFile', entry.name, 'transformType', transformType, ...
    'appliedOn', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'appliedBy', mfilename);

MotionCorrectionLog = rec([]);
if isfile(logFile)
    try
        S = load(logFile, 'MotionCorrectionLog');
        MotionCorrectionLog = S.MotionCorrectionLog;
    catch
        warning('Umitoolbox:correctMotionArtifact:UnreadableLog', ...
            'Could not read "%s"; it was replaced.', logFile);
    end
end
% One entry per file: a rewritten file supersedes its earlier entry.
MotionCorrectionLog(strcmpi({MotionCorrectionLog.file}, rec.file)) = [];
MotionCorrectionLog(end+1) = rec;
save(logFile, 'MotionCorrectionLog');

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
function [outData, corrAfter] = iApplyTransforms(entry, tforms, destPath)
%IAPPLYTRANSFORMS Warp every frame with its transform.
%
% Array input returns the corrected single array (DESTPATH is ignored). A
% .dat input is written to DESTPATH through a scratch file in the same
% folder (so the output only appears once the run has completed, and DESTPATH
% may be the input file itself) and its full path is returned. Also returns
% the correlation of each corrected frame with frame 1.

nt = entry.nt;
isDat = strcmp(entry.kind, 'dat');
outputView = imref2d([entry.ny, entry.nx]);
corrAfter = ones(nt, 1);

if isDat
    [destFolder, outStem, outExt] = fileparts(destPath);
    tmpPath = fullfile(destFolder, [outStem '_writing' outExt]);
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
function info = localPipelineInfo()
%LOCALPIPELINEINFO Return PipelineManager metadata for correctMotionArtifact.
%
% The step is a pipeline SETUP step: it has no DATA input or output and
% works on the folder (the batch form), so PipelineManager runs it before the
% processing steps. With Overwrite it corrects the imported channel files in
% place, and every later step reads the corrected files under the names
% ImportedChannels lists.

info = PipelineManager.createPipelineInfo( ...
    'correctMotionArtifact', ...
    ['Estimate per-frame motion against the first frame of every imported ' ...
     '.dat channel in the folder and correct it, optionally in place.']);

info = PipelineManager.addInput( ...
    info, ...
    'SaveFolder', ...
    'SaveFolder', ...
    ['Folder containing AcqInfos.mat and the imported .dat files; the ' ...
     'corrected files and the optional run files are written here.'], ...
    'kind', 'input', ...
    'position', 1, ...
    'callType', 'positional', ...
    'isData', false);

info = PipelineManager.addInput( ...
    info, ...
    'TransformType', ...
    'parameter', ...
    'Motion model: translation (fast) or similarity (rotation, scale and translation).', ...
    'kind', 'parameter', ...
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
    'callType', 'namevalue', ...
    'default', 100, ...
    'allowed', [1 Inf], ...
    'dataType', 'double');

info = PipelineManager.addInput( ...
    info, ...
    'Overwrite', ...
    'parameter', ...
    ['If true, every channel file is replaced in place by its corrected ' ...
     'version, so later steps use it; this cannot be undone. If false, each ' ...
     'is written as <name>_MotionCorrected.dat and later steps do not use it.'], ...
    'kind', 'parameter', ...
    'callType', 'namevalue', ...
    'default', false, ...
    'allowed', [true false], ...
    'dataType', 'logical');

info = PipelineManager.addInput( ...
    info, ...
    'AllowRecorrection', ...
    'parameter', ...
    ['If false, a channel that was already motion corrected is skipped (an ' ...
     'error if none is left). If true, it is corrected again after a warning.'], ...
    'kind', 'parameter', ...
    'callType', 'namevalue', ...
    'default', false, ...
    'allowed', [true false], ...
    'dataType', 'logical');

info = PipelineManager.addInput( ...
    info, ...
    'SaveShifts', ...
    'parameter', ...
    ['If true, saves each channel''s shift matrix (Nt-by-4 [tx, ty, rotationDeg, ' ...
     'scale]), correlations and run metadata to <name>_MotionCorrection.mat ' ...
     'and the QC figure to <name>_MotionCorrectionQC.png; nothing is ' ...
     'saved when false.'], ...
    'kind', 'parameter', ...
    'callType', 'namevalue', ...
    'default', true, ...
    'allowed', [true false], ...
    'dataType', 'logical');

info = PipelineManager.addInput( ...
    info, ...
    'ShowPlot', ...
    'parameter', ...
    'If true, shows the correlation and motion QC figure of each channel.', ...
    'kind', 'parameter', ...
    'callType', 'namevalue', ...
    'default', true, ...
    'allowed', [true false], ...
    'dataType', 'logical');

info = PipelineManager.addOutput( ...
    info, ...
    'correctedDatFiles', ...
    {'ImageTimeSeries','ProcessedData'}, ...
    'file', ...
    ['The imported channel .dat files: replaced in place when Overwrite is ' ...
     'true, otherwise joined by their <name>_MotionCorrected.dat copies. ' ...
     'MotionCorrectionLog.mat is written too; it is not declared as a ' ...
     'separate output because a non-returning file output cannot be optional.'], ...
    '*.dat', ...
    1, ...
    'isData', false, ...
    'returnsValue', false);

info.notes = { ...
    ['Uses the third-party DFTREGISTRATION algorithm (Guizar-Sicairos, ' ...
     'Thurman & Fienup, Opt. Lett. 33, 156-158, 2008), bundled as a ' ...
     'private helper.']; ...
    'Each channel is corrected independently against its own first frame.'; ...
    ['With Overwrite the uncorrected recordings are replaced; keep the raw ' ...
     'data to import again.']; ...
    ['A channel that was already motion corrected is skipped, and the step ' ...
     'fails when no channel is left to correct.']};

end
