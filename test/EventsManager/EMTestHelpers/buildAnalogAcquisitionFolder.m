function folderPath = buildAnalogAcquisitionFolder(folderPath, varargin)
%BUILDANALOGACQUISITIONFOLDER Create a synthetic analog acquisition folder.
%
% The folder contains a compatible info.txt and ai_00000.bin file with the
% trigger on StimAna1 and a valid camera trigger on channel 1.
%
% Name-Value Pairs:
%   'AISampleRate'   - Analog sample rate. Default = 10000.
%   'FrameRateHz'    - Imaging frame rate. Default = 60.
%   'NSamples'       - Number of samples before block padding. Default = 25000.
%   'PulseStarts'    - StimAna1 pulse starts in samples. Default = [5000 10000 15000].
%   'PulseWidth'     - Pulse width in samples. Default = 80.
%   'Amplitude'      - Pulse amplitude. Default = 5.
%   'CameraOn'       - Camera trigger start sample. Default = 1001.
%   'CameraOff'      - Camera trigger end sample. Default = 24000.
%   'StimName'       - Stimulation1 name. Default = 'Main'.
%
% Output:
%   folderPath - Folder containing synthetic acquisition files.

    p = inputParser;
    addRequired(p, 'folderPath', @(x) ischar(x) || isStringScalar(x));
    addParameter(p, 'AISampleRate', 10000, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'FrameRateHz', 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'NSamples', 25000, @(x) isnumeric(x) && isscalar(x) && x >= 10000);
    addParameter(p, 'PulseStarts', [5000 10000 15000], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'PulseWidth', 80, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Amplitude', 5, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'CameraOn', 1001, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'CameraOff', 24000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'StimName', 'Main', @(x) ischar(x) || isStringScalar(x));
    parse(p, folderPath, varargin{:});

    folderPath = convertStringsToChars(folderPath);
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end

    nSamples = p.Results.NSamples;
    nChan = 12;
    data = zeros(nSamples, nChan, 'single');

    cameraOn = p.Results.CameraOn;
    cameraOff = min(nSamples, p.Results.CameraOff);
    data(cameraOn:cameraOff, 1) = 5;

    stim = makePulseSignal(nSamples, p.Results.PulseStarts, p.Results.PulseWidth, ...
        'Amplitude', p.Results.Amplitude, 'Baseline', 0);
    data(:, 2) = stim;

    writeInfoTxtStimAna1(folderPath, ...
        'FrameRateHz', p.Results.FrameRateHz, ...
        'AISampleRate', p.Results.AISampleRate, ...
        'StimName', p.Results.StimName, ...
        'StimNRepeat', numel(p.Results.PulseStarts));
    writeSyntheticAIFile(folderPath, data, 'FileName', 'ai_00000.bin');
end
