function saveData(filename, data, varargin)
%SAVEDATA Save raw image-series data or derived UMIT data to disk.
%
%   saveData(filename, data)
%   saveData(filename, data, AcqInfoStream)
%   saveData(filename, data, AcqInfoStream, 'Append', true)
%
%   Inputs:
%       filename        - Full path of the file to be saved. The extension
%                         is forced based on the input data type.
%       data            - Either:
%                         1) non-empty single numeric 3D array, saved as .dat
%                         2) valid derived-data structure, saved as .umt
%
%   Optional inputs for numeric data only:
%       AcqInfoStream   - Acquisition info structure. Required only if
%                         AcqInfos.mat does not already exist in SaveFolder.
%
%   Name-Value options:
%       Append          - Logical scalar. If true, append numeric data to an
%                         existing .dat file. Default: false
%
%   Notes:
%       - Numeric 3D arrays are saved with ".dat".
%       - Structured derived data are saved as MAT-files with ".umt".
%       - If "AcqInfos.mat" is missing and numeric data are being saved,
%         AcqInfoStream must be provided.
%       - Numeric data must match the frame dimensions stored in
%         "AcqInfos.mat".

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'saveData';

addRequired(p, 'filename', @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(p, 'data', @(x) (isnumeric(x) && isa(x, 'single') && ~isempty(x)) || isstruct(x));
addOptional(p, 'AcqInfoStream', struct.empty(0,1), @isstruct);
addParameter(p, 'Append', false, @(x) islogical(x) && isscalar(x));

parse(p, filename, data, varargin{:});

filename = convertStringsToChars(filename);
[saveFolder, fileBase, ~] = fileparts(filename);

if isempty(saveFolder)
    saveFolder = pwd;
end

if isstruct(data)
    outFile = fullfile(saveFolder, [fileBase, '.umt']);
    save2umt(outFile, data);
else
    outFile = fullfile(saveFolder, [fileBase, '.dat']);
    save2dat(outFile, data, p.Results.AcqInfoStream, p.Results.Append);
end

disp(['Data saved as "' outFile '"']);
end

% =========================================================================
% Local functions
% =========================================================================

function save2dat(filePath, data, AcqInfoStream, b_append)
%SAVE2DAT Save numeric data array to a binary .dat file.

if ~(isnumeric(data) && isa(data, 'single') && ~isempty(data))
    error('Umitoolbox:saveData:invalidInput', ...
        'Numeric data must be a non-empty single array.');
end

if ndims(data) ~= 3
    error('Umitoolbox:saveData:invalidInput', ...
        'Numeric data must be a 3D single array.');
end

[saveFolder, fileName, ext] = fileparts(filePath);
if isempty(saveFolder)
    saveFolder = pwd;
end
if ~strcmpi(ext, '.dat')
    filePath = fullfile(saveFolder, [fileName, '.dat']);
end

if ~isfolder(saveFolder)
    error('Umitoolbox:saveData:invalidFolder', ...
        'Target folder does not exist: "%s".', saveFolder);
end

acqInfoFile = fullfile(saveFolder, 'AcqInfos.mat');

if ~isfile(acqInfoFile)
    if isempty(AcqInfoStream)
        error('Umitoolbox:saveData:missingAcqInfos', ...
            'Failed to save "%s". No "AcqInfos.mat" file found in folder "%s".', ...
            filePath, saveFolder);
    end

    save(acqInfoFile, 'AcqInfoStream');
end

S = load(acqInfoFile, 'AcqInfoStream');
if ~isfield(S, 'AcqInfoStream')
    error('Umitoolbox:saveData:invalidAcqInfos', ...
        '"AcqInfos.mat" does not contain variable "AcqInfoStream".');
end

AcqInfoStream = S.AcqInfoStream;

if ~isstruct(AcqInfoStream) || ~isscalar(AcqInfoStream) || ...
        ~isfield(AcqInfoStream, 'Height') || ~isfield(AcqInfoStream, 'Width')
    error('Umitoolbox:saveData:invalidAcqInfos', ...
        '"AcqInfoStream" must be a scalar struct containing "Height" and "Width".');
end

if ~isequal([size(data, 1), size(data, 2)], [AcqInfoStream.Height, AcqInfoStream.Width])
    error('Umitoolbox:saveData:invalidInput', ...
        ['Operation aborted. Data does not have the same frame dimensions ' ...
         'as in "AcqInfos.mat".']);
end

permission = 'w';
if b_append
    permission = 'a';
end

disp('Writing data to .DAT file ...');

fid = fopen(filePath, permission);
if fid == -1
    error('Umitoolbox:saveData:fileOpenFailed', ...
        'Could not open file for writing: "%s".', filePath);
end

cleanupObj = onCleanup(@() safeFclose(fid)); %#ok<NASGU>

nExpected = numel(data);
nWritten = fwrite(fid, data, 'single');

if nWritten ~= nExpected
    error('Umitoolbox:saveData:fileWriteFailed', ...
        'Failed to write all data to file "%s" (%d/%d elements written).', ...
        filePath, nWritten, nExpected);
end
end

function save2umt(filePath, data)
%SAVE2UMT Save derived-data structure to a .umt MAT-file.

validateUMTStruct(data);

disp('Writing data to .UMT file ...');
save(filePath, '-struct', 'data', '-mat');
end