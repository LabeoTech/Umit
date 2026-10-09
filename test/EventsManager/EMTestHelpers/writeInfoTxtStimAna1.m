function filePath = writeInfoTxtStimAna1(folderPath, varargin)
%WRITEINFOTXTSTIMANA1 Write a synthetic analog-stimulation info.txt file.
%
% The format is based on a real StimAna1 metadata example and is intended
% to be compatible with ReadInfoFile.
%
% Name-Value Pairs:
%   'FrameRateHz'      - Imaging frame rate. Default = 60.
%   'AISampleRate'     - Analog input sample rate. Default = 10000.
%   'AINChannels'      - Number of analog channels. Default = 12.
%   'Width'            - Image width. Default = 112.
%   'Height'           - Image height. Default = 112.
%   'StimName'         - Stimulation1 name. Default = 'Main'.
%   'StimPeriodMs'     - Stimulation1 period. Default = 333.
%   'StimDurationMs'   - Stimulation1 duration. Default = 5.
%   'StimPulseWidthMs' - Stimulation1 pulse width. Default = 167.
%   'StimNRepeat'      - Stimulation1 N Repeat. Default = 3.
%
% Output:
%   filePath - Full path to info.txt.

    p = inputParser;
    addRequired(p, 'folderPath', @(x) ischar(x) || isStringScalar(x));
    addParameter(p, 'FrameRateHz', 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'AISampleRate', 10000, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'AINChannels', 12, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Width', 112, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Height', 112, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'StimName', 'Main', @(x) ischar(x) || isStringScalar(x));
    addParameter(p, 'StimPeriodMs', 333, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'StimDurationMs', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'StimPulseWidthMs', 167, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'StimNRepeat', 3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(p, folderPath, varargin{:});

    folderPath = convertStringsToChars(folderPath);
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end

    filePath = fullfile(folderPath, 'info.txt');
    fid = fopen(filePath, 'w');
    assert(fid ~= -1, 'Failed to open %s for writing.', filePath);
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'Scan info\n');
    fprintf(fid, 'DateTime: 20230421_142936\n');
    fprintf(fid, 'Version: 4.3.4   2D Optogen plugin version: 3.3\n');
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
    fprintf(fid, 'AINChannels: %d\n', p.Results.AINChannels);
    chanNames = {'CameraTrig','StimAna1','StimAna2','AI1','AI2','AI3','AI4','AI5','AI6','AI7','AI8','CameraTrig2'};
    for ii = 1:p.Results.AINChannels
        fprintf(fid, 'AICh%d: %s\n', ii, chanNames{ii});
    end
    fprintf(fid, 'Illumination1: Red\n');
    fprintf(fid, 'Illumination2: Amber\n');
    fprintf(fid, 'Illumination3: Green\n');
    fprintf(fid, 'Stimulation: 1\n');
    fprintf(fid, 'Stimulation1 Name: %s\n', convertStringsToChars(p.Results.StimName));
    fprintf(fid, 'Stimulation1 Period: %g\n', p.Results.StimPeriodMs);
    fprintf(fid, 'Stimulation1 Duration: %g\n', p.Results.StimDurationMs);
    fprintf(fid, 'Stimulation1 Pulse Width: %g\n', p.Results.StimPulseWidthMs);
    fprintf(fid, 'Stimulation1 Burst Repeat: 1\n');
    fprintf(fid, 'Stimulation1 Burst Delay: 0\n');
    fprintf(fid, 'Stimulation1 Amplitude: 10\n');
    fprintf(fid, 'Stimulation1 Biphasic: 0\n');
    fprintf(fid, 'Stimulation1 Inter Stim: 10\n');
    fprintf(fid, 'Stimulation1 N Repeat: %d\n', p.Results.StimNRepeat);
    fprintf(fid, 'Stimulation1 Jitter: 3\n');
    fprintf(fid, 'Stimulation1 Output Channel: 1\n');
    fprintf(fid, 'Stimulation1 Test Run: 0\n');
    fprintf(fid, 'Stimulation2: 0\n');
end
