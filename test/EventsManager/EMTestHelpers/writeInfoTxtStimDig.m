function filePath = writeInfoTxtStimDig(folderPath, eventOrder, eventNames, eventDurationsSec, varargin)
%WRITEINFOTXTSTIMDIG Write a synthetic digital-stimulation info.txt file.
%
% The format is based on a real StimDig metadata example and is intended to
% be compatible with ReadInfoFile.
%
% Inputs:
%   folderPath         - Destination folder.
%   eventOrder         - Vector of stimulus IDs in presentation order.
%   eventNames         - Cell array of names indexed by stimulus ID.
%   eventDurationsSec  - Vector of durations in seconds indexed by stimulus ID.
%
% Name-Value Pairs:
%   'FrameRateHz'  - Imaging frame rate. Default = 120.
%   'AISampleRate' - Analog input sample rate. Default = 10000.
%   'Width'        - Image width. Default = 48.
%   'Height'       - Image height. Default = 48.
%
% Output:
%   filePath - Full path to info.txt.

    p = inputParser;
    addRequired(p, 'folderPath', @(x) ischar(x) || isStringScalar(x));
    addRequired(p, 'eventOrder', @(x) isnumeric(x) && isvector(x) && ~isempty(x));
    addRequired(p, 'eventNames', @(x) iscell(x) || isstring(x));
    addRequired(p, 'eventDurationsSec', @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'FrameRateHz', 120, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'AISampleRate', 10000, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Width', 48, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Height', 48, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, folderPath, eventOrder, eventNames, eventDurationsSec, varargin{:});

    eventOrder = p.Results.eventOrder(:)';
    eventNames = cellstr(p.Results.eventNames(:));
    eventDurationsSec = p.Results.eventDurationsSec(:);

    nIDs = max(eventOrder);
    assert(numel(eventNames) >= nIDs, 'eventNames must cover all event IDs.');
    assert(numel(eventDurationsSec) >= nIDs, 'eventDurationsSec must cover all event IDs.');

    folderPath = convertStringsToChars(folderPath);
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end

    filePath = fullfile(folderPath, 'info.txt');
    fid = fopen(filePath, 'w');
    assert(fid ~= -1, 'Failed to open %s for writing.', filePath);
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'Scan info\n');
    fprintf(fid, 'DateTime: 20230421_132932\n');
    fprintf(fid, 'Version: 4.3.3   2D Optogen plugin version: 3.3\n');
    fprintf(fid, 'Camera Model: CS2100M\n');
    fprintf(fid, 'FrameRateHz: %.6f\n', p.Results.FrameRateHz);
    fprintf(fid, 'Width: %d\n', p.Results.Width);
    fprintf(fid, 'Height: %d\n', p.Results.Height);
    fprintf(fid, 'Binning: 8\n');
    fprintf(fid, 'ExposureMsec: 0.833300\n');
    fprintf(fid, 'ExposureFluoMsec: 0.833300\n');
    fprintf(fid, 'ExposureSpeckleMsec: 0.833300\n');
    fprintf(fid, 'Rotation: 0\n');
    fprintf(fid, 'Delay first stim: 5.000000\n');
    fprintf(fid, 'AISampleRate: %d\n', p.Results.AISampleRate);
    fprintf(fid, 'AINChannels: 10\n');
    chanNames = {'CameraTrig','StimDig','AI1','AI2','AI3','AI4','AI5','AI6','AI7','AI8'};
    for ii = 1:numel(chanNames)
        fprintf(fid, 'AICh%d: %s\n', ii, chanNames{ii});
    end
    fprintf(fid, 'Illumination1: Red\n');
    fprintf(fid, 'Illumination2: Amber\n');
    fprintf(fid, 'Illumination3: Green\n');
    fprintf(fid, 'Stimulation: 2\n');
    fprintf(fid, 'Stimulation Repeat: %d\n', numel(eventOrder));
    fprintf(fid, 'Stimulation Randomize: Latin Square Design\n');
    fprintf(fid, 'Events Order:\t');
    fprintf(fid, '%d ', eventOrder);
    fprintf(fid, '\n\n');
    fprintf(fid, 'Events Description:\n');
    fprintf(fid, 'ID\tName\tCode\tDuration\tOpto Stim\tOpto Stim ID\n');
    for ii = 1:nIDs
        fprintf(fid, '%d\t%s\t%d\t%g\t0\t0\n', ii, eventNames{ii}, 100 + ii, eventDurationsSec(ii));
    end
end
