function run_ConvertToTiff(data,SaveFolder,varargin)
% RUN_CONVERTTOTIFF calls the function CONVERTTOTIFF from the IOI library (LabeoTech).
% In brief, this function creates a .TIFF file that can be opened in other
% softwares such as ImageJ.
% This function will create a .TIFF file for each imaging time series. In
% cases where the input data is a time series split by events, a .TIFF
% file will be created for each trial.
% Inputs:
%   data: numerical matrix containing imaging data.
%   SaveFolder: folder where the .TIFF file(s) will be saved.
%   Optional:
%   ByEvent (bool | default = FALSE): If TRUE, split the output
%   data by events. An "events.mat" file is necessary in the SaveFolder.

% Defaults:
default_opts = struct('ByEvent', false);
opts_values = struct('ByEvent', [true,false]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.

%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p,'data',@(x) isnumeric(x)); % Validate if the input is numerical matrix:
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,data,SaveFolder);
data = p.Results.data;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
clear p
% Create .TIF file name:
[~,outFileName,~] = fileparts(metaData.datFile);
if isempty(outFileName)
    outFileName = 'img_out';
end
% If the data is separated by Events, create one TIFF file per trial
if opts.ByEvent
    [data,condIndx,repIndx] = reshape_data_by_event(data,SaveFolder);          
    for ii = 1:size(data,1)
        str = [outFileName, '_', num2str(condIndx(ii)), 'c', num2str(repIndx(ii)),'r.tif'];
        ConvertToTiff(SaveFolder, squeeze(data(ii,:,:,:)), str)
    end
else
    ConvertToTiff(SaveFolder,data,[outFileName '.tif']);
end

end