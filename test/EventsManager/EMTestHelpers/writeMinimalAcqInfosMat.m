function acqInfo = writeMinimalAcqInfosMat(folderPath, varargin)
%WRITEMINIMALACQINFOSMAT Create a minimal AcqInfos.mat file for tests.
%
% Syntax:
%   acqInfo = writeMinimalAcqInfosMat(folderPath)
%   acqInfo = writeMinimalAcqInfosMat(folderPath, 'Name', Value, ...)
%
% Input:
%   folderPath - Destination folder.
%
% Name-Value Pairs:
%   'FrameRateHz'   - Imaging frame rate. Default = 60.
%   'AISampleRate'  - Analog input sample rate. Default = 10000.
%   'AINChannels'   - Number of analog channels. Default = 12.
%   'AIChanList'    - Cell array of analog channel names.
%   'StimName'      - Name for stimulation 1. Default = 'Main'.
%   'StimPeriodMs'  - Stimulation 1 period in ms. Default = 333.
%   'StimDurationMs'- Stimulation 1 duration in ms. Default = 5.
%
% Output:
%   acqInfo - Structure saved as AcqInfoStream inside AcqInfos.mat.

    p = inputParser;
    addRequired(p, 'folderPath', @(x) ischar(x) || isStringScalar(x));
    addParameter(p, 'FrameRateHz', 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'AISampleRate', 10000, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'AINChannels', 12, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'AIChanList', {}, @(x) iscell(x) || isstring(x));
    addParameter(p, 'StimName', 'Main', @(x) ischar(x) || isStringScalar(x));
    addParameter(p, 'StimPeriodMs', 333, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'StimDurationMs', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, folderPath, varargin{:});

    folderPath = convertStringsToChars(folderPath);
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end

    nChan = p.Results.AINChannels;
    aiChanList = p.Results.AIChanList;
    if isempty(aiChanList)
        defaultNames = {'CameraTrig','StimAna1','StimAna2','AI1','AI2','AI3','AI4','AI5','AI6','AI7','AI8','CameraTrig2'};
        assert(nChan <= numel(defaultNames), 'Provide AIChanList when AINChannels exceeds %d.', numel(defaultNames));
        aiChanList = defaultNames(1:nChan);
    else
        aiChanList = cellstr(aiChanList);
        assert(numel(aiChanList) == nChan, 'AIChanList must match AINChannels.');
    end

    acqInfo = struct();
    acqInfo.FrameRateHz = single(p.Results.FrameRateHz);
    acqInfo.AISampleRate = single(p.Results.AISampleRate);
    acqInfo.AINChannels = nChan;
    acqInfo.Width = 112;
    acqInfo.Height = 112;
    acqInfo.DateTime = '20260101_120000';
    for ii = 1:nChan
        acqInfo.(sprintf('AICh%d', ii)) = aiChanList{ii};
    end

    acqInfo.Stimulation1_Name = convertStringsToChars(p.Results.StimName);
    acqInfo.Stimulation1_Period = single(p.Results.StimPeriodMs);
    acqInfo.Stimulation1_Duration = single(p.Results.StimDurationMs);
    acqInfo.Stimulation1_Burst_Delay = single(0);

    AcqInfoStream = acqInfo; %#ok<NASGU>
    save(fullfile(folderPath, 'AcqInfos.mat'), 'AcqInfoStream');
end
