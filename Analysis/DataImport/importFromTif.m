function outFile = importFromTif(RawFolder,SaveFolder,varargin)
% IMPORTFROMTIF Imports single-channel imaging data stored in TIF format.
%
% This function reads TIFF files located in the specified "RawFolder" and
% their associated metadata from "info.json". It concatenates multi-file
% image sequences, applies optional spatial and temporal binning, and
% saves the processed data as .dat files in "SaveFolder".
%
% Parameters:
%   RawFolder   - Path to the folder containing the raw TIFF files.
%   SaveFolder  - Path to the folder where processed data will be saved.
%   opts        - (Optional) Struct with processing options:
%                 - 'BinningSpatial' : Spatial downsampling factor (default = 1)
%                 - 'BinningTemp'    : Temporal downsampling factor (default = 1)
%
% Returns:
%   outFile - Cell array containing the names of the saved .dat files.

% Defaults:
default_Output = {'fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'};%#ok This is here only as a reference for PIPELINEMANAGER.m. The real outputs will be stored in OUTFILE.
default_opts = struct('BinningSpatial', 1, 'BinningTemp', 1);
opts_values = struct('BinningSpatial', 2.^[0:5], 'BinningTemp',1:8);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
% Arguments validation:
p = inputParser;
addRequired(p, 'RawFolder', @isfolder);
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
parse(p, RawFolder, SaveFolder, varargin{:});
% Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
outFile = {};
clear p
%%%%

% Read info.json file in RawFolder:
assert(isfile(fullfile(RawFolder,'info.json')),'Failed to import TIF file(s). The "info.json" file is missing in "%s"!',RawFolder)
tif_metadata = jsondecode(fileread(fullfile(RawFolder,'info.json')));
if isstruct(tif_metadata.Tiffiles)
    % Force the structure Tiffiles to cell array for conformity. Doing this because if the
    % json contains different channels (TIF files) with different fielnames
    % it will be decoded as a cell array.
    tif_metadata.Tiffiles = arrayfun(@(x) x, tif_metadata.Tiffiles,'UniformOutput',false);
end
if isempty(tif_metadata.Tiffiles)
    error('Missing information about TIF files! Please check info.json file and try again!')
end

for ii = 1:length(tif_metadata.Tiffiles)
    % Look for multiple TIFs.
    % Here, we try to find Image sequences stored in multiple files with
    % the file name as prefix followed by a number. All files will be
    % concatenated into a single .dat file.
    tif_list = getTIFlist(RawFolder,tif_metadata.Tiffiles{ii}.filename);
    if isempty(tif_list)
        error('Data Import failed! No TIF found with name "%s" in folder "%s"!',tif_metadata.Tiffiles{ii}.filename,RawFolder)
    end
    % Get information from all files and preallocate the data:
    tif_info = cellfun(@imfinfo,tif_list, 'UniformOutput',false);
    tmp = imread(tif_info{1}(1).Filename,'Info',tif_info{1}(1));
    tot_data_sz = [size(tmp), sum(cellfun(@length,tif_info))];
    data = zeros(tot_data_sz,class(tmp));
    w = waitbar(0, 'Importing frames from TIFF files...');
    w.Children(1).Title.Interpreter = 'none';
    cnt = 1;
    for jj = 1:length(tif_list)
        this_tif_info = tif_info{jj};
        % Import frames:
        w.Children.Title.String = sprintf('Importing frames from file "%s"...',tif_list{jj});
        for kk = 1:length(this_tif_info)
            data(:,:,cnt) = imread(this_tif_info(1).Filename, 'Info', this_tif_info(kk));
            waitbar(kk/length(this_tif_info),w);
            cnt = cnt + 1;
        end
    end
    % Force data to "single" format :
    data = single(data);
    % Perform spatial and/or temporal binning:
    % Temporal Binning
    if( opts.BinningTemp > 1 )
        data = imresize3(data, [size(data,1), size(data,2),...
            size(data,3)/opts.BinningTemp], 'linear');
        % Update Temporal frequency:
        tif_metadata.Tiffiles{ii}.FrameRateHz = tif_metadata.Tiffiles{ii}.FrameRateHz/opts.BinningTemp;
    end
    % Spatial Binning
    if( opts.BinningSpatial > 1 )
        data = imresize(data,1/opts.BinningSpatial);
    end
    % Save everything to a .dat file in SaveFolder:
    datFileName = [lower(tif_metadata.Tiffiles{ii}.IlluminationColor) '.dat'];
    % Generate meta data:
    extraParams = tif_metadata.Tiffiles{ii};
    extraParams.tExposure = tif_metadata.Tiffiles{ii}.ExposureMsec;
    extraParams.Freq = tif_metadata.Tiffiles{ii}.FrameRateHz;
    fieldsToDelete = {'filename','FrameRateHz','ExposureMsec','IlluminationColor'};
    extraParams = rmfield(extraParams,fieldsToDelete);
    metaData = genMetaData(data,{'Y','X','T'}, extraParams);
    % Update datFile:
    metaData.datFile = fullfile(SaveFolder, datFileName);
    % Save data to .dat file:
    save2Dat(metaData.datFile, data, metaData);
    % Update "outFile" list:
    outFile = [outFile, {datFileName}];%#ok
    close(w);
end
% Create AcqInfo.mat file from global fields of the json file:
AcqInfoStream = rmfield(tif_metadata,'Tiffiles');
save(fullfile(SaveFolder,'AcqInfo.mat'), 'AcqInfoStream')
disp('Finished importFromTif.')
end
% Local functions
function tifNames = getTIFlist(folder,filename)
% GETTIFLIST Retrieves the list of TIFF files forming an image sequence.
%
% This function searches for multiple TIFF files that follow a common
% naming pattern (string followed by a number with or without a separator ("_" or "-"),
% indicating they are part of a sequential dataset.
%
% Parameters:
%   folder   - Path to the folder containing the TIFF files.
%   filename - Name of the primary TIFF file, used to infer a sequence.
%
% Returns:
%   tifNames - Cell array of ordered TIFF filenames that form the sequence.
%              Returns an empty array if no matching files are found.
%
% Errors:
%   - If no valid TIFF sequence is found, returns an empty cell array.

filename = convertStringsToChars(filename);
if ~endsWith(filename, '.tif','IgnoreCase',true)
    filename = [filename '.tif'];
end
prefix = regexp(filename, '^(.*?)([-_]\d+)?\.tif$','tokens','once','ignorecase');
prefix = prefix{1};

% Look for data stored in multiple files.
tifList = dir(fullfile(folder,[prefix '*.tif']));
if length(tifList) == 1  && strcmp(tifList.name,filename)
    % A single file was found and it is the given filename;
    tifNames = {tifList.name};
    return
elseif length(tifList) > 1
    % Multiple files were found:
    % Get list of files using the following rules:
    %  1) the prefix is followed by a common separator:
    %       "-","_", or none.
    %  2) a numeric sequence follows the prefix.
    
    % Get sequence:
    tifNames = {tifList.name};
    fprintf('Found sequence of files (N = %d) with prefix "%s"\nProcessing...\n',length(tifNames),prefix);
    exp = ['^' prefix '[-_]?(\d+)\.tif$'];
    sufixes = regexp(tifNames, exp, 'tokens','ignorecase');
    idx = ~cellfun(@isempty,sufixes);
    
    if ~any(idx) && any(strcmp({tifList.name},filename))
        % There were other similar files in the folder but they don't seem to
        % be part of a sequence. In this case return the filename.
        tifNames = filename;
        return
    end
    % Set file list in ascending numerical order
    tifNames = tifNames(idx);
    sufixes = sufixes(idx);sufixes = [sufixes{:}]; sufixes = [sufixes{:}];
    [~,fileOrd] = sort(cellfun(@str2double,sufixes));
    tifNames = tifNames(fileOrd);
else
    % No files found!
    tifNames = {};
end

end