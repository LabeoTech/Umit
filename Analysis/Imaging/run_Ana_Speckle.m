function [outData, metaData]= run_Ana_Speckle(SaveFolder, data, varargin)
% RUN_ANA_SPECKLE calls the function ANA_SPECKLE from the IOI library (LabeoTech).
% In brief, this function calculates blood flow (in arbitrary units) from 
% a Laser Speckle Contrast Imaging data

% Defaults:
default_Output = 'Flow.dat'; %#ok This line is here just for Pipeline management
default_opts = struct('bNormalize', false, 'SpeckleFileName','speckle');
opts_values = struct('bNormalize',[false, true], 'SpeckleFileName',{{'speckle'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
addRequired(p,'data',@(x) ( isnumeric(x) && ndims(x)==3 ) || ischar(x));
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,SaveFolder,data, varargin{:});
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
clear p
% Run Ana_Speckle function from IOI library:
disp('Calculating  blood flow...')
if ischar(data)
    b_RAMsafeMode = true;
else
    b_RAMsafeMode = false;
end
[outData, metaData] = Ana_Speckle(SaveFolder,opts.bNormalize, opts.SpeckleFileName,'bRAMsafe',b_RAMsafeMode);
if b_RAMsafeMode
    % Update meta data with extra params from original file
    originalMetaData = load(fullfile(SaveFolder,strrep(data,'.dat','.mat')));
    fnToMerge = setdiff(fieldnames(originalMetaData),fieldnames(metaData));
    for ii = 1:length(fnToMerge)
        metaData.(fnToMerge{ii}) = originalMetaData.(fnToMerge{ii});
    end
end
    
disp('Finished Speckle Mapping.')
end
