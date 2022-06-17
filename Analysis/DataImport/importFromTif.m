function outFile = importFromTif(RawFolder,SaveFolder,varargin)
% IMPORTFROMTIF is a data import function that reads single-channel imaging
% data stored in TIF format. The function looks for any .tif (e.g. fluo.tif)
% file inside the "RawFolder" and the associate info file (e.g. fluo.txt)
% containing the recording meta data such as sample rate and exposure time.

% Defaults:
default_Output = {'fluo.dat', 'red.dat', 'green.dat', 'yellow.dat', 'speckle.dat'};%#ok This is here only as a reference for PIPELINEMANAGER.m. The real outputs will be stored in OUTFILE.
default_opts = struct('BinningSpatial', 1, 'BinningTemp', 1);
opts_values = struct('BinningSpatial', 2.^[0:6], 'BinningTemp',2.^[0:6]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
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

% Read all TIF and TXT files inside "RawFolder":
tif_list = [dir(fullfile(RawFolder,'*.tif'));dir(fullfile(RawFolder,'*.tiff'))];
txt_list = dir(fullfile(RawFolder,'*.txt'));
[~,tifNames,~] = cellfun(@fileparts,{tif_list.name}, 'UniformOutput',false);
[~,txtNames,~] = cellfun(@fileparts,{txt_list.name}, 'UniformOutput',false);
% Select those that have a .txt file with the same name: 
[idx,locB] = ismember(tifNames,txtNames);
if ~any(idx) || isempty(idx)
    error('umIToolbox:importFromTif:FileNotFound','The RawFolder does not contain any valid .tif files!');
elseif ~all(idx)
    fprintf('The following .tif files do not have the associate .txt files and will be ignored:\n');
    fprintf('%s\n',tifNames{~idx})
    fprintf('-----------\n');
end      
tif_list = tif_list(idx);
txt_list = txt_list(locB(locB~=0));
% Read TIF files and each TXT file:
disp('Importing data from TIFF files...')
for i = 1:length(tif_list)
    tif_file = fullfile(tif_list(i).folder, tif_list(i).name);
    info_tif = imfinfo(tif_file);        
    % Check if data has a valid format:
    tmp = imread(tif_file, 'Info', info(1));
    if ~ismember(class(tmp),{'uint8', 'uint16', 'uint32'})
        error('umIToolbox:importFromTif:InvalidFile',...
            ['Invalid data format: This function accepts only unsigned integer data '...
            'with 8, 16 or 32 Bits.'])
    end
    % Import frames:
    disp(['[Started] Data import from file: ' tif_list(i).name])
    data = zeros([size(tmp),length(info_tif)],class(tmp));
    for j = 1:length(info_tif)
        data(:,:,j) = imread(tif_file, 'Info', info(j));
    end
    disp(['[Completed] Data import from file: ' tif_list(i).name])
    % Read TXT file:
    disp(['[Started] Reading meta data from file: ' tif_list(i).name])
    
    %%% TO BE CONTINUED..%%%
    
    % Perform spatial and/or temporal binning:
    
    %%% TO BE CONTINUED..%%%
    
    % Save everything to a .dat file in SaveFolder:
    
    %%% TO BE CONTINUED..%%%
    
end
    
end
