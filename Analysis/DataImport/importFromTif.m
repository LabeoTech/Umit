function outFile = importFromTif(RawFolder, SaveFolder, varargin)
% IMPORTFROMTIF Import single-channel imaging data stored in TIFF format.
%
% This function reads TIFF files located in "RawFolder" and their
% associated metadata from "info.json". It concatenates multi-file image
% sequences, applies optional spatial and temporal binning, and saves the
% processed data as .dat files in "SaveFolder".
%
% Syntax:
%   outFile = importFromTif(RawFolder, SaveFolder)
%   outFile = importFromTif(RawFolder, SaveFolder, 'Name', Value, ...)
%   info    = importFromTif('pipelineInfo')
%
% Inputs:
%   RawFolder   - Path to the folder containing the raw TIFF files.
%   SaveFolder  - Path to the folder where processed .dat files will be
%                 saved.
%
% Name-Value parameters:
%   BinningSpatial - Spatial downsampling factor.
%                    Allowed values: 1, 2, 4, 8, 16, 32
%   BinningTemp    - Temporal downsampling factor.
%                    Allowed values: 1:8
%
% Output:
%   outFile     - Cell array containing the names of the saved .dat files.
%
% Notes:
%   - This function creates data outputs as file-backed .dat files, so its
%     pipelineInfo output is declared with outputMode='file' and isData=true.
%   - Folder-level metadata are stored in "AcqInfos.mat" using GENACQINFO.
%   - All imported channels within the same SaveFolder must have identical
%     dimensions and acquisition timing metadata.

% -------------------------------------------------------------------------
% Shared declarations
% -------------------------------------------------------------------------
defaultOutput = {'fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'};
allowedBinningSpatial = 2.^[0:5];
allowedBinningTemp = 1:8;

% -------------------------------------------------------------------------
% pipelineInfo query
% -------------------------------------------------------------------------
if nargin == 1 && (ischar(RawFolder) || (isstring(RawFolder) && isscalar(RawFolder))) ...
        && strcmpi(strtrim(char(string(RawFolder))), 'pipelineInfo')
    outFile = localGetPipelineInfo();
    return
end

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'RawFolder', @isfolder);
addRequired(p, 'SaveFolder', @isfolder);
addParameter(p, 'BinningSpatial', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && any(x == allowedBinningSpatial));
addParameter(p, 'BinningTemp', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && any(x == allowedBinningTemp));
parse(p, RawFolder, SaveFolder, varargin{:});

RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
BinningSpatial = p.Results.BinningSpatial;
BinningTemp = p.Results.BinningTemp;

clear p

outFile = {};

% -------------------------------------------------------------------------
% Read and validate info.json
% -------------------------------------------------------------------------
infoPath = fullfile(RawFolder, 'info.json');
assert(isfile(infoPath), ...
    'Failed to import TIFF file(s). The "info.json" file is missing in "%s"!', RawFolder);

tif_metadata = jsondecode(fileread(infoPath));

assert(isfield(tif_metadata, 'Tiffiles'), ...
    'Failed to import TIFF file(s). The "info.json" file does not contain the "Tiffiles" field.');

if isstruct(tif_metadata.Tiffiles)
    tif_metadata.Tiffiles = num2cell(tif_metadata.Tiffiles);
end

assert(iscell(tif_metadata.Tiffiles) && ~isempty(tif_metadata.Tiffiles), ...
    'Missing information about TIFF files. Please check info.json and try again.');

% -------------------------------------------------------------------------
% Prepare folder-level AcqInfoStream
% -------------------------------------------------------------------------
acqInfoPath = fullfile(SaveFolder, 'AcqInfos.mat');
AcqInfoStream = [];

if isfile(acqInfoPath)
    tmp = load(acqInfoPath);

    if isfield(tmp, 'AcqInfoStream')
        AcqInfoStream = tmp.AcqInfoStream;
    elseif isfield(tmp, 'AcqInfos')
        AcqInfoStream = tmp.AcqInfos;
    else
        fn = fieldnames(tmp);
        assert(~isempty(fn), ...
            'The existing "AcqInfos.mat" file in "%s" does not contain any variables.', SaveFolder);
        AcqInfoStream = tmp.(fn{1});
    end

    assert(isstruct(AcqInfoStream) && isscalar(AcqInfoStream), ...
        'The existing "AcqInfos.mat" file in "%s" does not contain a valid acquisition-info structure.', SaveFolder);
end

extraInfo = rmfield(tif_metadata, 'Tiffiles');
coreFields = intersect(fieldnames(extraInfo), {'Width','Height','Length','FrameRateHz','ExposureMsec'});
if ~isempty(coreFields)
    extraInfo = rmfield(extraInfo, coreFields);
end
extraInfo.BinningSpatial = BinningSpatial;
extraInfo.BinningTemp = BinningTemp;

% -------------------------------------------------------------------------
% Import each TIFF channel/sequence
% -------------------------------------------------------------------------
for ii = 1:numel(tif_metadata.Tiffiles)
    tifEntry = tif_metadata.Tiffiles{ii};

    assert(isstruct(tifEntry), ...
        'Invalid TIFF entry found in "info.json". Each Tiffiles entry must be a structure.');
    assert(isfield(tifEntry, 'filename') && ~isempty(tifEntry.filename), ...
        'Invalid TIFF entry found in "info.json". Missing field "filename".');
    assert(isfield(tifEntry, 'FrameRateHz') && ~isempty(tifEntry.FrameRateHz), ...
        'Invalid TIFF entry found in "info.json". Missing field "FrameRateHz" for "%s".', ...
        char(string(tifEntry.filename)));
    assert(isfield(tifEntry, 'ExposureMsec') && ~isempty(tifEntry.ExposureMsec), ...
        'Invalid TIFF entry found in "info.json". Missing field "ExposureMsec" for "%s".', ...
        char(string(tifEntry.filename)));
    assert(isfield(tifEntry, 'IlluminationColor') && ~isempty(tifEntry.IlluminationColor), ...
        'Invalid TIFF entry found in "info.json". Missing field "IlluminationColor" for "%s".', ...
        char(string(tifEntry.filename)));

    % Resolve TIFF sequence
    tif_list = getTIFlist(RawFolder, tifEntry.filename);
    if isempty(tif_list)
        error('Data import failed. No TIFF file found with name "%s" in folder "%s".', ...
            char(string(tifEntry.filename)), RawFolder);
    end

    % Read TIFF metadata from all files in the sequence
    tif_info = cellfun(@(x) imfinfo(fullfile(RawFolder, x)), tif_list, 'UniformOutput', false);

    % Read first frame to determine data class and dimensions
    firstFrame = imread(tif_info{1}(1).Filename, 'Info', tif_info{1}(1));
    totalFrames = sum(cellfun(@numel, tif_info));
    data = zeros([size(firstFrame,1), size(firstFrame,2), totalFrames], class(firstFrame));

    % Import full sequence
    w = waitbar(0, 'Importing frames from TIFF files...', 'Name', 'importFromTif');
    try
        cnt = 1;
        for jj = 1:numel(tif_list)
            this_tif_info = tif_info{jj};

            for kk = 1:numel(this_tif_info)
                data(:,:,cnt) = imread(this_tif_info(1).Filename, 'Info', this_tif_info(kk));

                waitbar( ...
                    ((jj - 1) + (kk / numel(this_tif_info))) / numel(tif_list), ...
                    w, ...
                    sprintf('Importing frames from file "%s"...', tif_list{jj}));

                cnt = cnt + 1;
            end
        end

        if isgraphics(w)
            close(w);
        end
    catch ME
        if isgraphics(w)
            close(w);
        end
        rethrow(ME)
    end

    % Force imported data to single precision
    data = single(data);

    % Temporal binning
    frameRateHz = tifEntry.FrameRateHz;
    if BinningTemp > 1
        assert(mod(size(data,3), BinningTemp) == 0, ...
            ['Temporal binning failed for "%s". The number of frames (%d) must be ' ...
             'divisible by BinningTemp (%d).'], ...
            char(string(tifEntry.filename)), size(data,3), BinningTemp);

        newLength = size(data,3) / BinningTemp;
        data = imresize3(data, [size(data,1), size(data,2), newLength], 'linear');
        frameRateHz = frameRateHz / BinningTemp;
    end

    % Spatial binning
    if BinningSpatial > 1
        data = imresize(data, 1 / BinningSpatial);
    end

    data = single(data);

    % Build candidate AcqInfoStream for this imported file
    thisAcqInfo = genAcqInfo( ...
        size(data,2), ...
        size(data,1), ...
        size(data,3), ...
        frameRateHz, ...
        tifEntry.ExposureMsec, ...
        extraInfo);

    if isempty(AcqInfoStream)
        AcqInfoStream = thisAcqInfo;
    else
        assert(isfield(AcqInfoStream, 'Width') && isfield(AcqInfoStream, 'Height') && ...
               isfield(AcqInfoStream, 'Length') && isfield(AcqInfoStream, 'FrameRateHz') && ...
               isfield(AcqInfoStream, 'ExposureMsec'), ...
            ['The existing "AcqInfos.mat" file in "%s" does not contain the required ' ...
             'fields Width, Height, Length, FrameRateHz, and ExposureMsec.'], SaveFolder);

        assert(isequal(AcqInfoStream.Width, thisAcqInfo.Width) && ...
               isequal(AcqInfoStream.Height, thisAcqInfo.Height) && ...
               isequal(AcqInfoStream.Length, thisAcqInfo.Length) && ...
               isequal(AcqInfoStream.FrameRateHz, thisAcqInfo.FrameRateHz) && ...
               isequal(AcqInfoStream.ExposureMsec, thisAcqInfo.ExposureMsec), ...
            ['Imported TIFF files in "%s" do not share the same dimensions and acquisition timing. ' ...
             'All generated .dat files within a SaveFolder must have identical Width, Height, Length, ' ...
             'FrameRateHz, and ExposureMsec.'], SaveFolder);
    end

    % Save data to .dat file
    datFileName = [lower(char(string(tifEntry.IlluminationColor))) '.dat'];
    datFilePath = fullfile(SaveFolder, datFileName);

    fid = fopen(datFilePath, 'w');
    assert(fid ~= -1, 'Failed to create data file "%s".', datFilePath);
    c = onCleanup(@() fclose(fid));

    nWritten = fwrite(fid, data, 'single');
    assert(nWritten == numel(data), ...
        'Failed to write all data to "%s". Expected %d elements, wrote %d.', ...
        datFilePath, numel(data), nWritten);

    clear c

    outFile = [outFile, {datFileName}]; %#ok<AGROW>
end

% -------------------------------------------------------------------------
% Save folder-level acquisition metadata
% -------------------------------------------------------------------------
save(acqInfoPath, 'AcqInfoStream');

disp('Finished importFromTif.')

% -------------------------------------------------------------------------
% Local pipelineInfo factory
% -------------------------------------------------------------------------
function info = localGetPipelineInfo()

    info = PipelineManager.createPipelineInfo( ...
        mfilename, ...
        'Import single-channel imaging data stored in TIFF format.');

    info = PipelineManager.addInput( ...
        info, ...
        'RawFolder', ...
        'RawFolder', ...
        'Path to the folder containing the raw TIFF files.', ...
        'kind', 'input', ...
        'position', 1, ...
        'callType', 'positional', ...
        'isData', false);

    info = PipelineManager.addInput( ...
        info, ...
        'SaveFolder', ...
        'SaveFolder', ...
        'Path to the folder where processed .dat files will be saved.', ...
        'kind', 'input', ...
        'position', 2, ...
        'callType', 'positional', ...
        'isData', false);

    info = PipelineManager.addInput( ...
        info, ...
        'BinningSpatial', ...
        'parameter', ...
        'Spatial downsampling factor.', ...
        'kind', 'parameter', ...
        'default', 1, ...
        'allowed', allowedBinningSpatial, ...
        'callType', 'namevalue');

    info = PipelineManager.addInput( ...
        info, ...
        'BinningTemp', ...
        'parameter', ...
        'Temporal downsampling factor.', ...
        'kind', 'parameter', ...
        'default', 1, ...
        'allowed', allowedBinningTemp, ...
        'callType', 'namevalue');

    info = PipelineManager.addOutput( ...
        info, ...
        'outFile', ...
        'Unknown', ...
        'file', ...
        'Imported TIFF data file(s) saved in SaveFolder.', ...
        defaultOutput, ...
        1, ...
        'isData', true, ...
        'saveFileName', '');
end

end

% -------------------------------------------------------------------------
% Local function
% -------------------------------------------------------------------------
function tifNames = getTIFlist(folder, filename)
% GETTIFLIST Retrieve the ordered list of TIFF files forming an image sequence.
%
% This function searches for TIFF files that follow a common naming pattern
% where the base filename may be followed by a numeric suffix with or
% without a separator ("_" or "-").
%
% Inputs:
%   folder   - Path to the folder containing the TIFF files.
%   filename - Name of the primary TIFF file used to infer the sequence.
%
% Output:
%   tifNames - Cell array of ordered TIFF filenames that form the sequence.
%              Returns {} if no matching files are found.

filename = convertStringsToChars(filename);

if ~endsWith(filename, '.tif', 'IgnoreCase', true)
    filename = [filename '.tif'];
end

tok = regexp(filename, '^(.*?)([-_]\d+)?\.tif$', 'tokens', 'once', 'ignorecase');
assert(~isempty(tok), ...
    'Invalid TIFF filename "%s". Expected a .tif file name.', filename);

prefix = tok{1};
prefixEsc = regexptranslate('escape', prefix);

% Look for data stored in multiple files
tifList = dir(fullfile(folder, [prefix '*.tif']));

if numel(tifList) == 1 && strcmpi(tifList.name, filename)
    tifNames = {tifList.name};
    return
elseif numel(tifList) > 1
    tifNames = {tifList.name};
    fprintf('Found sequence of files (N = %d) with prefix "%s"\nProcessing...\n', ...
        numel(tifNames), prefix);

    expr = ['^' prefixEsc '[-_]?(\d+)\.tif$'];
    suffixes = regexp(tifNames, expr, 'tokens', 'ignorecase');
    idx = ~cellfun(@isempty, suffixes);

    if ~any(idx) && any(strcmpi({tifList.name}, filename))
        % Similar files exist but do not belong to a numeric sequence.
        tifNames = {filename};
        return
    end

    tifNames = tifNames(idx);
    suffixes = suffixes(idx);
    suffixes = [suffixes{:}];
    suffixes = [suffixes{:}];

    [~, fileOrd] = sort(cellfun(@str2double, suffixes));
    tifNames = tifNames(fileOrd);
else
    tifNames = {};
end

end