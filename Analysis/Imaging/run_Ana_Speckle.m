function [outData, metaData]= run_Ana_Speckle(SaveFolder, varargin)
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
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,SaveFolder, varargin{:});
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
clear p
% Run Ana_Speckle function from IOI library:
disp('Calculating  blood flow...')
[outData, metaData] = Ana_Speckle(SaveFolder,opts.bNormalize, opts.SpeckleFileName);
disp('Finished Speckle Mapping.')
end
