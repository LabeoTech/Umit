function [tf, reason, marker] = isUmitDataset(folder)
%ISUMITDATASET Format-neutral test: does FOLDER hold imported imaging data?
%
%   [tf, reason, marker] = isUmitDataset(folder)
%
%   A folder is a UMIT dataset when it holds at least one of:
%     1. a .dat file with a valid self-describing header (readDatHeader);
%     2. an AcqInfos.mat with a scalar AcqInfoStream struct (the import
%        record of folders imported before the .dat header, and of
%        imports that write it today).
%   Neither marker is required on its own, so the test keeps working when
%   AcqInfos.mat moves back to raw-data scope (.dat header paradigm shift).
%   Read-only.
%
%   Outputs:
%       tf      - logical scalar
%       reason  - '' when tf is true, else why the folder is not a dataset
%       marker  - 'datHeader', 'acqInfos' or ''
%
%   Used by UMITProjectStore (addSession and session validation) instead
%   of requiring AcqInfos.mat (audit section 12.1 item 5).
%
%   See also readDatHeader, isDatWithHeader.

folder = char(string(folder));
tf = false;
marker = '';
if ~isfolder(folder)
    reason = sprintf('the folder does not exist: %s', folder);
    return
end

datFiles = dir(fullfile(folder, '*.dat'));
datFiles = datFiles(~[datFiles.isdir]);
for iFile = 1:numel(datFiles)
    datPath = fullfile(folder, datFiles(iFile).name);
    try
        if isDatWithHeader(datPath)
            readDatHeader(datPath);
            tf = true;
            reason = '';
            marker = 'datHeader';
            return
        end
    catch
        % An unreadable or corrupt header is not a marker; keep looking.
    end
end

acqInfoPath = fullfile(folder, 'AcqInfos.mat');
if ~isfile(acqInfoPath)
    reason = sprintf(['no imported data were found in "%s": it has no .dat file ' ...
        'with a UMIT header and no AcqInfos.mat'], folder);
    return
end
try
    loaded = struct();
    if ismember('AcqInfoStream', who('-file', acqInfoPath))
        loaded = load(acqInfoPath, 'AcqInfoStream', '-mat');
    end
catch ME
    reason = sprintf('AcqInfos.mat in "%s" could not be read: %s', folder, ME.message);
    return
end
if ~isfield(loaded, 'AcqInfoStream') || ~isstruct(loaded.AcqInfoStream) || ...
        ~isscalar(loaded.AcqInfoStream)
    reason = sprintf(['AcqInfos.mat in "%s" must contain a scalar AcqInfoStream ' ...
        'struct, and the folder has no .dat file with a UMIT header'], folder);
    return
end
tf = true;
reason = '';
marker = 'acqInfos';
end
