function run_ConvertToTiff(data, metaData, SaveFolder)
% RUN_CONVERTTOTIFF calls the function
% CONVERTTOTIFF from the IOI library (LabeoTech).
% In brief, this function creates a .TIFF file that can be opened in other
% softwares such as ImageJ.
% This function will create a .TIFF file for each imaging time series. In
% cases where the input data is a time series split by events, a .TIFF
% file will be created for each trial.
% Inputs:
%   data: numerical matrix containing imaging data.
%   metaData: .mat file with meta data associated with "data".
%   SaveFolder: folder where the .TIFF file(s) will be saved.
%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p,'data',@(x) isnumeric(x) | ischar(x)); % Validate if the input is numerical matrix:
addRequired(p,'metaData', @(x) isa(x,'matlab.io.MatFile') | isstruct(x)); % MetaData associated to "data".
addRequired(p, 'SaveFolder', @isfolder);
% Parse inputs:
parse(p,data, metaData, SaveFolder);
data = p.Results.data;
metaData = p.Results.metaData;
SaveFolder = p.Results.SaveFolder;
clear p
% Create .TIF file name:
[~,outFileName,~] = fileparts(metaData.datFile);
if isempty(outFileName)
    outFileName = 'img_out';
else
    outFileName = ['img_' outFileName];
end

if ischar(data)
    [~,filename,ext] = fileparts(data);
    filename = [filename,ext];
    b_fromFile = true;
    
else
    b_fromFile = false;
end

% If the data is separated by Events, create one TIFF file per trial
if any(strcmpi(metaData.dim_names, 'E'))
    if b_fromFile
        mapFile = mapDatFile(fullfile(SaveFolder,filename));
    end
    % Create list of suffix to identify conditions and trials:
    IDlist = unique(metaData.eventID);
    str = cell(size(metaData.eventID));
    for i = 1:length(IDlist)
        indx = find(metaData.eventID == IDlist(i));
        tmp = arrayfun(@(x) [outFileName '_C' num2str(i) '_R' num2str(x) '.tif'], indx,'UniformOutput',false);
        str(indx) = tmp;
    end
    disp('Creating one TIF file per trial...');
    for i = 1:size(data,1)
        if b_fromFile
            tmp = squeeze(mapFile.Data.data(i,:,:,:));
            ConvertToTiff(SaveFolder, tmp, str{i});
        else
            ConvertToTiff(SaveFolder, squeeze(data(i,:,:,:)), str{i})
        end
    end
else
    if b_fromFile
        ConvertToTiff(SaveFolder, filename);
    else
        ConvertToTiff(SaveFolder, data, [outFileName '.tif']);
    end
end

end