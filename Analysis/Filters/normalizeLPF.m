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
%       1) Numeric array with dimensions Y x X x T
%       2) Numeric array with dimensions Y x X x T x E (event-split data)
%       3) Filename to a .dat file with axes Y-X-T or Y-X-T-E
%       UMT structs and .umt/.mat files are not supported.
%
%   Input/output behavior:
%       - If the input is a numeric array, the output is a numeric array
%         with the same size.
%       - If the input is a .dat filename, the output is a .dat filename
%         ("normLPF.dat" in SaveFolder) with the same axes and sizes.
%
%   Inputs:
%       data       - Input data in one of the accepted forms above.
%       SaveFolder - Output folder (for .dat inputs).
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
%       FrameRateHz      - Frame rate of DATA (Hz), used for the filter.
%                          PipelineManager injects it from the data; a
%                          .dat input's header provides it otherwise.
%                          In-RAM arrays need it explicitly; AcqInfos.mat
%                          is not used.
%
%   Output:
%       outData     - Filtered data with the same representation type as
%                     the input.
%
%   Notes:
%       - The filtering algorithm itself is delegated to the IOI library
%         function "NormalisationFiltering". This wrapper only handles
%         input resolution, validation, and low-RAM orchestration.
%       - Raw .dat input is passed through to NormalisationFiltering in its
%         own file mode; the chunking itself happens inside
%         NormalisationFiltering, not in a low-RAM helper in this wrapper.
%         Y-X-T and Y-X-T-E files are streamed in X slabs, so the recording
%         is never loaded whole.
%       - Event-split data (an E axis) are filtered trial-by-trial along T
%         and keep the E dimension unchanged; with bApplyExpFit each trial
%         gets its own exponential fit.
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
addParameter(p, 'FrameRateHz', []);

parse(p, data, SaveFolder, varargin{:});

SaveFolder = char(string(p.Results.SaveFolder));
BaselineCutoffHz = double(p.Results.BaselineCutoffHz);
SignalCutoffHz = double(p.Results.SignalCutoffHz);
bNormalize = p.Results.Normalize;
bApplyExpFit = p.Results.bApplyExpFit;
explicitRate = p.Results.FrameRateHz;

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

% Frame rate of the data itself: the explicit FrameRateHz (injected by
% PipelineManager), else the .dat header. AcqInfos.mat is not used
% (resolveDataInfoValue).
Fs = [];
if isFileInput && strcmp(ext, '.dat')
    % The layout is checked first: an axis layout without time (for example
    % Y-X-E) has no frame rate to resolve.
    assertDatLayout(loadMetaData(dataFile), ...
        {{'Y','X','T'}, {'Y','X','T','E'}}, 'normalizeLPF');
    Fs = resolveDataInfoValue('frameRateHz', explicitRate, dataFile, mfilename);
elseif isnumeric(data) || islogical(data)
    Fs = resolveDataInfoValue('frameRateHz', explicitRate, data, mfilename);
end
if ~isempty(Fs)
    iCheckCutoffs(Fs, BaselineCutoffHz, SignalCutoffHz);
end

% -------------------------------------------------------------------------
% Case 1: YXT or YXTE array in RAM
% -------------------------------------------------------------------------
if isnumeric(data) || islogical(data)

    validateattributes(data, {'numeric','logical'}, {'nonempty'}, ...
        mfilename, 'data');
    if ~(ndims(data) == 3 || ndims(data) == 4)
        error('normalizeLPF:InvalidArrayInput', ...
            'Numeric input must be YXT or YXTE.');
    end

    if ndims(data) == 3
        outData = iFilterArray(single(data), BaselineCutoffHz, SignalCutoffHz, ...
            bNormalize, bApplyExpFit, Fs);
    else
        % Event-split data: every trial is filtered along T on its own.
        outData = zeros(size(data), 'single');
        for iTrial = 1:size(data, 4)
            outData(:,:,:,iTrial) = iFilterArray(single(data(:,:,:,iTrial)), ...
                BaselineCutoffHz, SignalCutoffHz, bNormalize, bApplyExpFit, Fs);
        end
    end
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

        otherwise
            error('normalizeLPF:UnsupportedInputFile', ...
                'Unsupported input file extension "%s". Only .dat files are supported.', ext);
    end
end

error('normalizeLPF:UnsupportedInputType', ...
    ['Input "data" must be a YXT or YXTE array or a .dat filename. ' ...
     'UMT structs and .umt/.mat files are not supported.']);

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
            ['Input data. Accepted forms: YXT or YXTE array, or a .dat ' ...
             'filename with axes Y-X-T or Y-X-T-E.'], ...
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
            'Output folder.', ...
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

        info = PipelineManager.addInput( ...
            info, ...
            'FrameRateHz', ...
            'sourceInfo', ...
            'Frame rate of the input data (Hz), injected from the data.', ...
            'kind', 'sourceInfo', ...
            'sourceField', 'frameRateHz', ...
            'required', false);

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
% Helper: cutoff checks
% =========================================================================
function iCheckCutoffs(Fs, BaselineCutoffHz, SignalCutoffHz)
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
end
