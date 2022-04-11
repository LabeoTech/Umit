function [outData, metaData]= run_Ana_Speckle(SaveFolder)
% RUN_ANA_SPECKLE calls the function ANA_SPECKLE from the IOI library (LabeoTech).
% In brief, this function calculates blood flow (in arbitrary units) from 
% a Laser Speckle Contrast Imaging data

% Defaults:
default_Output = 'Flow.dat'; % This line is here just for Pipeline management...
                               % In this case, the string must be inside a cell
                               % in order to PipelineManager to work.

%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Parse inputs:
parse(p,SaveFolder);
SaveFolder = p.Results.SaveFolder;
clear p
% Run Ana_Speckle function from IOI library:
disp('Calculating  blood flow...')
[outData, metaData] = Ana_Speckle(SaveFolder);
disp('Finished Speckle Mapping.')
end
