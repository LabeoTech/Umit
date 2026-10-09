function folderPath = buildDigitalAcquisitionFolder(folderPath, varargin)
%BUILDDIGITALACQUISITIONFOLDER Create a synthetic digital acquisition folder.
%
% The folder contains a compatible info.txt and ai_00000.bin file with the
% trigger on StimDig and a valid camera trigger on channel 1.
%
% Name-Value Pairs:
%   'AISampleRate' - Analog sample rate. Default = 10000.
%   'FrameRateHz'  - Imaging frame rate. Default = 120.
%   'NSamples'     - Number of samples before block padding. Default = 22000.
%   'EventOrder'   - Stimulus ID presentation order. Default = [3 1 2].
%   'EventNames'   - Cell array of event names indexed by ID.
%   'Durations'    - Durations in seconds indexed by ID.
%   'PulseStarts'  - Digital onset sample indices. Default = [5000 9000 13000].
%   'PulseWidth'   - Pulse width in samples. Default = 40.
%   'CameraOn'     - Camera trigger start sample. Default = 1001.
%   'CameraOff'    - Camera trigger end sample. Default = 21000.
%
% Output:
%   folderPath - Folder containing synthetic acquisition files.

    p = inputParser;
    addRequired(p, 'folderPath', @(x) ischar(x) || isStringScalar(x));
    addParameter(p, 'AISampleRate', 10000, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'FrameRateHz', 120, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'NSamples', 22000, @(x) isnumeric(x) && isscalar(x) && x >= 10000);
    addParameter(p, 'EventOrder', [3 1 2], @(x) isnumeric(x) && isvector(x) && ~isempty(x));
    addParameter(p, 'EventNames', {'block1','block2','block3'}, @(x) iscell(x) || isstring(x));
    addParameter(p, 'Durations', [3 5 1], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'PulseStarts', [5000 9000 13000], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'PulseWidth', 40, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'CameraOn', 1001, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'CameraOff', 21000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    parse(p, folderPath, varargin{:});

    folderPath = convertStringsToChars(folderPath);
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end

    nSamples = p.Results.NSamples;
    nChan = 10;
    data = zeros(nSamples, nChan, 'single');

    cameraOn = p.Results.CameraOn;
    cameraOff = min(nSamples, p.Results.CameraOff);
    data(cameraOn:cameraOff, 1) = 5;

    dig = makePulseSignal(nSamples, p.Results.PulseStarts, p.Results.PulseWidth, ...
        'Amplitude', 5, 'Baseline', 0);
    data(:, 2) = dig;

    writeInfoTxtStimDig(folderPath, p.Results.EventOrder, p.Results.EventNames, p.Results.Durations, ...
        'FrameRateHz', p.Results.FrameRateHz, 'AISampleRate', p.Results.AISampleRate);
    writeSyntheticAIFile(folderPath, data, 'FileName', 'ai_00000.bin');
end
