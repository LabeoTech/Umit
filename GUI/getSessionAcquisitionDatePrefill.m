function acquisitionDate = getSessionAcquisitionDatePrefill(saveFolder)
%GETSESSIONACQUISITIONDATEPREFILL Read an optional session date prefill.
%
%   acquisitionDate = getSessionAcquisitionDatePrefill(saveFolder)
%
% Reads only the canonical AcqInfos.mat file. A scalar
% AcqInfoStream.DateTime value in yyyyMMdd_HHmmss format is converted to
% a date-only datetime. Missing or malformed metadata returns NaT.

acquisitionDate = NaT;
metadataFile = fullfile(char(string(saveFolder)), 'AcqInfos.mat');
if ~isfile(metadataFile)
    return
end

try
    loaded = load(metadataFile, 'AcqInfoStream', '-mat');
    if ~isfield(loaded, 'AcqInfoStream') || ...
            ~isstruct(loaded.AcqInfoStream) || ...
            ~isscalar(loaded.AcqInfoStream) || ...
            ~isfield(loaded.AcqInfoStream, 'DateTime')
        return
    end

    dateText = string(loaded.AcqInfoStream.DateTime);
    if ~isscalar(dateText) || ...
            strlength(strtrim(dateText)) == 0
        return
    end

    parsedDate = datetime(strtrim(dateText), ...
        'InputFormat', 'yyyyMMdd_HHmmss');
    if ~isnat(parsedDate)
        acquisitionDate = dateshift(parsedDate, 'start', 'day');
    end
catch
    acquisitionDate = NaT;
end
end
