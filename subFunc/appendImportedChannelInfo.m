function AcqInfoStream = appendImportedChannelInfo(target, channelInfo, varargin)
%APPENDIMPORTEDCHANNELINFO Append imported-channel metadata to AcqInfos.mat.
%
%   AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, channelInfo)
%   AcqInfoStream = appendImportedChannelInfo(acqInfoPath, channelInfo)
%   AcqInfoStream = appendImportedChannelInfo(saveFolder, channelInfo)
%   AcqInfoStream = appendImportedChannelInfo(..., 'Overwrite', true)
%
%   Adds or updates AcqInfoStream.ImportedChannels entries for .dat files
%   created by import functions. This is not a folder file ledger. It only
%   records the image time-series files created during data import.
%
%   Required channelInfo fields:
%       DatFile     - Imported .dat file name.
%       Length      - Number of frames in the imported .dat file.
%       FrameRateHz - Frame rate of the imported .dat file.
%
%   Optional channelInfo fields:
%       Tag         - Short channel tag, e.g. 'red', 'fluo_475'.
%       Color       - Original/imported color label.
%       ExposureMsec- Exposure time in milliseconds.
%       CamIdx      - Camera index for multicamera imports. Default: 1 for
%                     non-multicamera acquisitions.
%
%   Notes:
%       - Spatial metadata remain at the top level of AcqInfoStream.
%       - RepeatCount and SequenceIdx are intentionally not stored.
%       - If target is a path, the updated AcqInfoStream is saved back to
%         AcqInfos.mat.

p = inputParser;
p.FunctionName = 'appendImportedChannelInfo';

addRequired(p, 'target', @(x) isstruct(x) || ischar(x) || isstring(x));
addRequired(p, 'channelInfo', @(x) isstruct(x) && ~isempty(x));
addParameter(p, 'Overwrite', false, @(x) islogical(x) && isscalar(x));

parse(p, target, channelInfo, varargin{:});

overwrite = p.Results.Overwrite;
acqInfoPath = '';

if isstruct(target)
    AcqInfoStream = target;
    if ~isscalar(AcqInfoStream)
        error('Umitoolbox:appendImportedChannelInfo:invalidInput', ...
            'AcqInfoStream input must be a scalar struct.');
    end
else
    targetPath = char(string(target));

    if isfolder(targetPath)
        acqInfoPath = fullfile(targetPath, 'AcqInfos.mat');
    else
        [folderPath, fileBase, ext] = fileparts(targetPath);
        if isempty(folderPath)
            folderPath = pwd;
        end
        if isempty(ext)
            acqInfoPath = fullfile(folderPath, fileBase, 'AcqInfos.mat');
        else
            acqInfoPath = fullfile(folderPath, [fileBase, ext]);
        end
    end

    if ~isfile(acqInfoPath)
        error('Umitoolbox:appendImportedChannelInfo:fileNotFound', ...
            'AcqInfos.mat file not found: "%s".', acqInfoPath);
    end

    S = load(acqInfoPath, 'AcqInfoStream');
    if ~isfield(S, 'AcqInfoStream')
        error('Umitoolbox:appendImportedChannelInfo:invalidAcqInfos', ...
            'File "%s" does not contain variable AcqInfoStream.', acqInfoPath);
    end
    AcqInfoStream = S.AcqInfoStream;
end

if ~isstruct(AcqInfoStream) || ~isscalar(AcqInfoStream)
    error('Umitoolbox:appendImportedChannelInfo:invalidAcqInfos', ...
        'AcqInfoStream must be a scalar struct.');
end

newEntries = iNormalizeChannelInfo(channelInfo, AcqInfoStream);

if isfield(AcqInfoStream, 'ImportedChannels') && ...
        ~isempty(AcqInfoStream.ImportedChannels)
    importedChannels = iNormalizeChannelInfo(AcqInfoStream.ImportedChannels, AcqInfoStream);
else
    importedChannels = iEmptyImportedChannels();
end

for iEntry = 1:numel(newEntries)
    thisEntry = newEntries(iEntry);
    idxExisting = find(strcmpi({importedChannels.DatFile}, thisEntry.DatFile), 1, 'first');

    if isempty(idxExisting)
        importedChannels(end+1) = thisEntry; %#ok<AGROW>
    elseif overwrite
        importedChannels(idxExisting) = thisEntry;
    elseif ~isequaln(importedChannels(idxExisting), thisEntry)
        error('Umitoolbox:appendImportedChannelInfo:duplicateChannel', ...
            ['Imported channel "%s" already exists with different metadata. ' ...
             'Use Overwrite=true to replace it.'], ...
            thisEntry.DatFile);
    end
end

AcqInfoStream.ImportedChannels = importedChannels;

if ~isempty(acqInfoPath)
    save(acqInfoPath, 'AcqInfoStream');
end

end

% =========================================================================
% Local helpers
% =========================================================================
function out = iNormalizeChannelInfo(channelInfo, AcqInfoStream)
%INORMALIZECHANNELINFO Normalize imported-channel metadata entries.

out = iEmptyImportedChannels();
channelInfo = channelInfo(:).';

for iEntry = 1:numel(channelInfo)
    raw = channelInfo(iEntry);

    requiredFields = {'DatFile', 'Length', 'FrameRateHz'};
    missingFields = setdiff(requiredFields, fieldnames(raw));
    if ~isempty(missingFields)
        error('Umitoolbox:appendImportedChannelInfo:missingField', ...
            'Imported channel metadata are missing field(s): %s.', ...
            strjoin(missingFields, ', '));
    end

    datFile = char(string(raw.DatFile));
    [~, datBase, datExt] = fileparts(datFile);
    if isempty(datExt)
        datExt = '.dat';
    end
    datFile = [datBase, datExt];

    if ~strcmpi(datExt, '.dat')
        error('Umitoolbox:appendImportedChannelInfo:invalidDatFile', ...
            'Imported channel DatFile must have extension .dat.');
    end

    len = double(raw.Length);
    freq = double(raw.FrameRateHz);

    validateattributes(len, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive', 'integer'}, ...
        'appendImportedChannelInfo', 'Length');
    validateattributes(freq, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, ...
        'appendImportedChannelInfo', 'FrameRateHz');

    if isfield(raw, 'Tag') && ~isempty(raw.Tag)
        tag = char(string(raw.Tag));
    elseif isfield(raw, 'Color') && ~isempty(raw.Color)
        tag = iNormalizeTag(raw.Color);
    else
        tag = iNormalizeTag(datBase);
    end

    if isfield(raw, 'Color') && ~isempty(raw.Color)
        color = char(string(raw.Color));
    else
        color = tag;
    end

    if isfield(raw, 'ExposureMsec') && ~isempty(raw.ExposureMsec)
        exposureMsec = double(raw.ExposureMsec);
    else
        exposureMsec = [];
    end

    if isfield(raw, 'CamIdx') && ~isempty(raw.CamIdx)
        camIdx = double(raw.CamIdx);
    elseif ~isfield(AcqInfoStream, 'MultiCam') || ~logical(AcqInfoStream.MultiCam)
        camIdx = 1;
    else
        camIdx = [];
    end

    out(end+1) = struct( ... %#ok<AGROW>
        'DatFile', datFile, ...
        'Tag', tag, ...
        'Color', color, ...
        'Length', len, ...
        'FrameRateHz', freq, ...
        'ExposureMsec', exposureMsec, ...
        'CamIdx', camIdx);
end

end

function out = iEmptyImportedChannels()
%IEMPTYIMPORTEDCHANNELS Return an empty ImportedChannels struct array.

out = struct( ...
    'DatFile', {}, ...
    'Tag', {}, ...
    'Color', {}, ...
    'Length', {}, ...
    'FrameRateHz', {}, ...
    'ExposureMsec', {}, ...
    'CamIdx', {});

end

function tag = iNormalizeTag(value)
%INORMALIZETAG Create a simple lowercase imported-channel tag.

tag = lower(char(string(value)));
tag = strrep(tag, ' ', '_');
tag = strrep(tag, '.dat', '');
if strcmpi(tag, 'amber')
    tag = 'yellow';
end

end
