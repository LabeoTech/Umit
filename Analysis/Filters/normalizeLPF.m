function outData = normalizeLPF(data, SaveFolder, varargin)
%NORMALIZELPF Normalize image data by low-pass filtering.
%
%   outData = normalizeLPF(data, SaveFolder)
%   outData = normalizeLPF(data, SaveFolder, ...
%       'BaselineCutoffHz', baselineCutoffHz, ...
%       'SignalCutoffHz', signalCutoffHz, ...
%       'Normalize', tfNormalize, ...
%       'bApplyExpFit', tfApplyExpFit)
%
%   This function wraps the IOI library function "NormalisationFiltering".
%   The filtering algorithm consists in creating two low-passed versions of
%   the signal with cut-off frequencies "BaselineCutoffHz" and
%   "SignalCutoffHz", then subtracting them. Optionally, the result can be
%   normalized by the baseline component to express the signal as DeltaR/R.
%
%   Accepted input forms:
%       1) Numeric 3-D array with dimensions Y x X x T
%       2) Filename to a .dat file storing Y x X x T data
%       3) UMT struct
%       4) Filename to a .umt or .mat file containing a UMT struct
%
%   Input/output behavior:
%       - If the input is a numeric 3-D array, the output is a numeric
%         3-D array with the same size.
%       - If the input is a .dat filename, the output is a .dat filename.
%       - If the input is a UMT struct, the output is a UMT struct.
%       - If the input is a .umt or .mat filename, the file is loaded in
%         RAM and the output is a UMT struct.
%
%   Inputs:
%       data       - Input data in one of the accepted forms above.
%       SaveFolder - Folder containing AcqInfos.mat used to retrieve the
%                    frame rate for filter validation and execution.
%
%   Name-Value parameters:
%       BaselineCutoffHz - Low cut-off frequency used to estimate the slow
%                          baseline component. Default: 0.0083
%
%       SignalCutoffHz   - Higher cut-off frequency used to preserve the
%                          signal component. Default: 1
%
%       Normalize        - Logical scalar. If true, express the filtered
%                          signal as DeltaR/R. Default: true
%
%       bApplyExpFit     - Logical scalar. If true, apply exponential decay
%                          correction inside NormalisationFiltering.
%                          Default: false
%
%   Output:
%       outData     - Filtered data with the same representation type as
%                     the input.
%
%   Notes:
%       - The filtering algorithm itself is delegated to the IOI library
%         function "NormalisationFiltering". This wrapper only handles
%         input resolution, validation, low-RAM orchestration, and UMT I/O.
%       - Raw .dat input is passed through to NormalisationFiltering in its
%         own file mode; the chunking itself happens inside
%         NormalisationFiltering, not in a low-RAM helper in this wrapper.
%       - UMT input must have kind = 'image'.
%       - All UMT entries must use dimensions:
%             {'Y','X','T'} or {'Y','X','T','E'}
%       - Entries with an E dimension are filtered trial-by-trial and keep
%         the E dimension unchanged.
%       - If a .umt or .mat file is provided, RAM-safe mode is not
%         available and the UMT content is loaded into RAM.
%       - For in-RAM (non-.dat) inputs, NaN pixels are replaced by 0 before
%         filtering and restored afterward. This biases the low-pass
%         baseline near mask borders, and with Normalize=true can feed a
%         near-zero divisor there.

default_Output = 'normLPF.dat'; %#ok<NASGU>

if nargin == 1 && (ischar(data) || (isstring(data) && isscalar(data))) ...
        && strcmpi(strtrim(char(string(data))), 'pipelineInfo')
    outData = localPipelineInfo();
    return
end

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'data');
addRequired(p, 'SaveFolder', @(x) ischar(x) || (isstring(x) && isscalar(x)));

addParameter(p, 'BaselineCutoffHz', 0.0083, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'SignalCutoffHz', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'Normalize', true, ...
    @(x) islogical(x) && isscalar(x));
addParameter(p, 'bApplyExpFit', false, ...
    @(x) islogical(x) && isscalar(x));

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
BaselineCutoffHz = double(p.Results.BaselineCutoffHz);
SignalCutoffHz = double(p.Results.SignalCutoffHz);
bNormalize = p.Results.Normalize;
bApplyExpFit = p.Results.bApplyExpFit;

if ~isfolder(SaveFolder)
    error('normalizeLPF:InvalidSaveFolder', ...
        'SaveFolder "%s" does not exist.', SaveFolder);
end

% Resolve a file input's path/extension up front (if any) so the frame
% rate below is read from the file actually being processed, not an
% arbitrary file in SaveFolder.
isFileInput = ischar(data) || (isstring(data) && isscalar(data));
dataFile = '';
ext = '';
if isFileInput
    dataFile = char(string(data));
    if ~isfile(dataFile)
        altPath = fullfile(SaveFolder, dataFile);
        if isfile(altPath)
            dataFile = altPath;
        else
            error('normalizeLPF:InputFileNotFound', ...
                'Input file "%s" was not found.', data);
        end
    end
    [~,~,ext] = fileparts(dataFile);
    ext = lower(ext);
end

if isFileInput && ismember(ext, {'.dat','.umt'})
    Fs = iGetFrameRateHz(SaveFolder, dataFile);
else
    Fs = iGetFrameRateHz(SaveFolder);
end

if BaselineCutoffHz < 0 || BaselineCutoffHz > Fs/2
    error('normalizeLPF:InvalidCutoff', ...
        'BaselineCutoffHz must be between 0 and the Nyquist frequency.');
end

if SignalCutoffHz <= 0 || SignalCutoffHz > Fs/2
    error('normalizeLPF:InvalidCutoff', ...
        'SignalCutoffHz must be > 0 and <= the Nyquist frequency.');
end

if SignalCutoffHz < BaselineCutoffHz
    error('normalizeLPF:InvalidCutoff', ...
        'SignalCutoffHz must be >= BaselineCutoffHz.');
end

% -------------------------------------------------------------------------
% Case 1: Raw YXT array in RAM
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)

    validateattributes(data, {'numeric','logical'}, {'nonempty','3d'}, ...
        mfilename, 'data');

    outData = iFilterArray(single(data), BaselineCutoffHz, SignalCutoffHz, ...
        bNormalize, bApplyExpFit, Fs);
    return
end

% -------------------------------------------------------------------------
% Case 2: File input
% -------------------------------------------------------------------------
if isFileInput

    switch ext
        case '.dat'
            outData = normalizeLPF_lowRAMmode( ...
                dataFile, SaveFolder, ...
                BaselineCutoffHz, SignalCutoffHz, ...
                bNormalize, bApplyExpFit, Fs);
            return

        case {'.umt','.mat'}
            warning('normalizeLPF:UMTFileLoadsInRAM', ...
                ['RAM-safe mode is not available for data stored in this format. ' ...
                 'Loading the UMT content into RAM.']);
            data = iLoadUMTFromFile(dataFile);

        otherwise
            error('normalizeLPF:UnsupportedInputFile', ...
                'Unsupported input file extension "%s".', ext);
    end
end

% -------------------------------------------------------------------------
% Case 3: UMT struct in RAM
% -------------------------------------------------------------------------
if ~isstruct(data)
    error('normalizeLPF:UnsupportedInputType', ...
        ['Input "data" must be a YXT array, a .dat filename, ' ...
         'a UMT struct, or a .umt/.mat filename containing a UMT struct.']);
end

[entryNames, entryData, entryDims, sourceLabels, sourceEventInfo, hasE] = ...
    iExtractValidUMTData(data);

out = data;
out.version = data.version;
out.kind = data.kind;
out.data = struct();

for iEntry = 1:numel(entryNames)

    value = entryData{iEntry};
    dimNames = entryDims{iEntry};

    if isequal(dimNames, {'Y','X','T'})
        filtData = iFilterArray(value, BaselineCutoffHz, SignalCutoffHz, ...
            bNormalize, bApplyExpFit, Fs);

    elseif isequal(dimNames, {'Y','X','T','E'})
        nTrials = size(value, 4);
        filtData = zeros(size(value), 'like', value);

        for iTrial = 1:nTrials
            trial = value(:,:,:,iTrial);
            filtData(:,:,:,iTrial) = iFilterArray(trial, ...
                BaselineCutoffHz, SignalCutoffHz, ...
                bNormalize, bApplyExpFit, Fs);
        end

    else
        error('normalizeLPF:InvalidUMTEntryDims', ...
            ['Entry "%s" must use dimNames {''Y'',''X'',''T''} or ' ...
             '{''Y'',''X'',''T'',''E''}.'], ...
            entryNames{iEntry});
    end

    out.data.(entryNames{iEntry}) = struct( ...
        'value', filtData, ...
        'dimNames', {dimNames});
end

if ~isempty(fieldnames(sourceLabels))
    out.labels = sourceLabels;
elseif isfield(out, 'labels')
    out = rmfield(out, 'labels');
end

if any(hasE)
    out = appendUMTEventInfo(out, ...
        'eventID', sourceEventInfo.eventID, ...
        'repetitionIndex', sourceEventInfo.repetitionIndex, ...
        'eventName', sourceEventInfo.eventName, ...
        'eventAxisMode', sourceEventInfo.eventAxisMode, ...
        'overwrite', true);
else
    if isfield(out, 'eventInfo')
        out = rmfield(out, 'eventInfo');
    end
    validateUMTStruct(out, 'requireEventInfo', true);
end

outData = out;

% =========================================================================
% Local pipeline info
% =========================================================================
    function info = localPipelineInfo()
        info = PipelineManager.createPipelineInfo(mfilename, ...
            ['Normalize image data by low-pass filtering using ' ...
             'NormalisationFiltering.']);

        info.version = '1.0.0';

        info = PipelineManager.addInput( ...
            info, ...
            'data', ...
            {'ImageTimeSeries','ProcessedData','UnknownDataType'}, ...
            ['Input data. Accepted forms: YXT array, .dat filename, ' ...
             'UMT struct, or .umt file containing one UMT struct.'], ...
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
            'Folder containing AcqInfos.mat.', ...
            'kind', 'input', ...
            'position', 2, ...
            'callType', 'positional', ...
            'isData', false);

        info = PipelineManager.addInput( ...
            info, ...
            'BaselineCutoffHz', ...
            'parameter', ...
            'Low cut-off frequency used to estimate the slow baseline component.', ...
            'kind', 'parameter', ...
            'default', 0.0083, ...
            'allowed', [0 Inf], ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'SignalCutoffHz', ...
            'parameter', ...
            'Higher cut-off frequency used to preserve the signal component.', ...
            'kind', 'parameter', ...
            'default', 1, ...
            'allowed', [0 Inf], ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'Normalize', ...
            'parameter', ...
            'If true, express the filtered signal as DeltaR/R.', ...
            'kind', 'parameter', ...
            'default', true, ...
            'allowed', [true false], ...
            'callType', 'namevalue');

        info = PipelineManager.addInput( ...
            info, ...
            'bApplyExpFit', ...
            'parameter', ...
            'If true, apply exponential decay correction.', ...
            'kind', 'parameter', ...
            'default', false, ...
            'allowed', [true false], ...
            'callType', 'namevalue');

        info = PipelineManager.addOutput( ...
            info, ...
            'outData', ...
            {'ImageTimeSeries','ProcessedData'}, ...
            'data', ...
            'Low-pass filtered output.', ...
            'normLPF.dat', ...
            1, ...
            'isData', true);
    end
end

% =========================================================================
% Helper: Low-RAM execution for raw .dat input
% =========================================================================
function outFile = normalizeLPF_lowRAMmode( ...
    inFile, SaveFolder, BaselineCutoffHz, SignalCutoffHz, ...
    bNormalize, bApplyExpFit, Fs)
%NORMALIZELPF_LOWRAMMODE Execute normalizeLPF on a raw .dat file.
%
% This helper preserves the original file-based behavior by delegating the
% filtering to NormalisationFiltering in file mode and returning the output
% filename.

outFile = fullfile(SaveFolder, 'normLPF.dat');

NormalisationFiltering( ...
    pwd, inFile, ...
    BaselineCutoffHz, ...
    SignalCutoffHz, ...
    bNormalize, ...
    bApplyExpFit, ...
    Fs, outFile);

end

% =========================================================================
% Helper: Filter one in-memory array while preserving NaNs
% =========================================================================
function outArray = iFilterArray(inArray, BaselineCutoffHz, SignalCutoffHz, ...
    bNormalize, bApplyExpFit, Fs)
%IFILTERARRAY Apply NormalisationFiltering to one numeric YXT block.

if ~isa(inArray, 'single')
    workArray = single(inArray);
else
    workArray = inArray;
end

idxNaN = isnan(workArray);
if any(idxNaN(:))
    workArray(idxNaN) = 0;
end

outArray = NormalisationFiltering( ...
    pwd, workArray, ...
    BaselineCutoffHz, ...
    SignalCutoffHz, ...
    bNormalize, ...
    bApplyExpFit, ...
    Fs);

outArray = single(outArray);

if any(idxNaN(:))
    outArray(idxNaN) = NaN;
end

end

% =========================================================================
% Helper: Extract and validate image-backed data from a UMT structure
% =========================================================================
function [entryNames, entryData, entryDims, labels, eventInfo, hasE] = iExtractValidUMTData(umt)

validateUMTStruct(umt, 'requireEventInfo', false);

if ~strcmpi(umt.kind, 'image')
    error('normalizeLPF:InvalidUMTKind', ...
        ['Operation aborted. UMT input must have kind = "image". ' ...
         'This function does not support non-image UMT structures.']);
end

entryNames = fieldnames(umt.data);
if isempty(entryNames)
    error('normalizeLPF:EmptyUMTData', ...
        'Operation aborted. UMT data is empty.');
end

entryData = cell(size(entryNames));
entryDims = cell(size(entryNames));
hasE = false(size(entryNames));

for iEntry = 1:numel(entryNames)
    thisEntry = umt.data.(entryNames{iEntry});
    thisDims = cellstr(string(thisEntry.dimNames));

    if ~(isequal(thisDims, {'Y','X','T'}) || isequal(thisDims, {'Y','X','T','E'}))
        error('normalizeLPF:InvalidUMTEntry', ...
            ['Operation aborted. All entries in the input UMT must use ' ...
             'dimNames {''Y'',''X'',''T''} or {''Y'',''X'',''T'',''E''}.' ...
             '\nInvalid entry: "%s".'], ...
            entryNames{iEntry});
    end

    entryData{iEntry} = thisEntry.value;
    entryDims{iEntry} = thisDims;
    hasE(iEntry) = isequal(thisDims, {'Y','X','T','E'});
end

if isfield(umt, 'labels')
    labels = umt.labels;
else
    labels = struct();
end

if any(hasE)
    if ~isfield(umt, 'eventInfo')
        error('normalizeLPF:MissingEventInfo', ...
            ['Operation aborted. The input UMT contains entries with an E ' ...
             'dimension but has no shared top-level eventInfo.']);
    end
    eventInfo = umt.eventInfo;
else
    eventInfo = struct();
end

end

% =========================================================================
% Helper: Load UMT from file
% =========================================================================
function umt = iLoadUMTFromFile(filePath)

[~,~,ext] = fileparts(filePath);
ext = lower(ext);

switch ext
    case '.umt'
        try
            tmp = loadData(filePath);
            if isstruct(tmp) && isscalar(tmp) && ...
                    all(ismember({'version','kind','data'}, fieldnames(tmp)))
                umt = tmp;
                return
            end
        catch
        end
        S = load(filePath, '-mat');

    case '.mat'
        S = load(filePath);

    otherwise
        error('normalizeLPF:InvalidUMTFile', ...
            'Unsupported UMT file extension "%s".', ext);
end

fn = fieldnames(S);
for iField = 1:numel(fn)
    candidate = S.(fn{iField});
    if isstruct(candidate) && isscalar(candidate) && ...
            all(ismember({'version','kind','data'}, fieldnames(candidate)))
        umt = candidate;
        return
    end
end

error('normalizeLPF:NoUMTFoundInFile', ...
    'No scalar UMT struct was found in "%s".', filePath);

end

% =========================================================================
% Helper: Get frame rate from AcqInfos.mat
% =========================================================================
function freqHz = iGetFrameRateHz(SaveFolder, dataFile)
%IGETFRAMERATEHZ Resolve the frame rate for the data being processed.
%
% When dataFile is given, resolve it from that file's own metadata via
% loadMetaData -- never from an arbitrary, unrelated file in SaveFolder.
% Otherwise (a raw numeric array or an in-RAM UMT struct, neither of which
% has a file identity of its own) fall back to the single authoritative
% AcqInfos.mat.

if nargin >= 2 && ~isempty(dataFile)
    meta = loadMetaData(dataFile);
    if ~isfield(meta, 'Freq') || isempty(meta.Freq)
        error('normalizeLPF:MissingFrameRate', ...
            'loadMetaData did not return Freq for "%s".', dataFile);
    end
    freqHz = double(meta.Freq);
    return
end

acqInfoFile = fullfile(SaveFolder, 'AcqInfos.mat');
if ~isfile(acqInfoFile)
    error('normalizeLPF:MissingReferenceData', ...
        'Could not determine frame rate because "AcqInfos.mat" was not found in "%s".', ...
        SaveFolder);
end

S = load(acqInfoFile, 'AcqInfoStream');
if ~isfield(S, 'AcqInfoStream') || ~isfield(S.AcqInfoStream, 'FrameRateHz') || ...
        isempty(S.AcqInfoStream.FrameRateHz)
    error('normalizeLPF:MissingFrameRate', ...
        '"AcqInfos.mat" in "%s" does not define FrameRateHz.', SaveFolder);
end

freqHz = double(S.AcqInfoStream.FrameRateHz);

end
