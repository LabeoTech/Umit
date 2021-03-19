function outFile = run_ImagesClassification(RawFolder, SaveFolder, varargin)
% RUN_IMAGESCLASSIFICATION calls the function
% IMAGESCLASSIFICATION from the IOI library (LabeoTech).

%%% Arguments parsing and validation %%%
p = inputParser;
% Raw folder:
addRequired(p, 'RawFolder', @isfolder);
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure:
default_opts = struct('BinningSpatial', 1, 'BinningTemp', 1, 'b_SubROI', false, 'b_IgnoreStim', true);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Output File:
default_Output = {'fChan_475.dat','fChan.dat', 'rChan.dat', 'gChan.dat', 'yChan.dat'}; % This is here only as a reference for PIPELINEMANAGER.m . The real outputs will be stored in OUTFILE.
addOptional(p, 'Output', default_Output, @(x) ischar(x) || isstring(x) || iscell(x));
parse(p, RawFolder, SaveFolder, varargin{:});
%Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
outFile = {};
%%%%
% Calls function from IOI library. Temporary for now.
ImagesClassification(RawFolder, opts.BinningSpatial, opts.BinningTemp, opts.b_SubROI, opts.b_IgnoreStim);
cd(RawFolder)
chanList = dir('*Chan*.dat'); chanList = {chanList.name};
for i = 1:length(chanList)
    chanName = chanList{i};
    switch chanName
        case 'rChan.dat'
            MetaDataFileName = 'Data_red.mat';
        case 'yChan.dat'
            MetaDataFileName = 'Data_yellow.mat';            
        case 'gChan.dat'
            MetaDataFileName = 'Data_green.mat';
        case {'fChan.dat', 'fChan_475.dat'}
            MetaDataFileName = 'Data_Fluo.mat';            
    end
    a =  matfile(fullfile(RawFolder,MetaDataFileName), 'Writable', true);
    a.fileUUID = char(java.util.UUID.randomUUID);
    a.Datatype = 'single';
    a.datName = 'data';
    % TEMPORARY FIX.
    filePath = fullfile(SaveFolder,chanName);
    newmDfile = strrep(filePath, '.dat', '_info.mat');
    statusDat = movefile(chanName, filePath);
    statusMat = movefile(MetaDataFileName, newmDfile);
    if ~statusDat
        disp(['Failed to transfer ' chanName ' to ' SaveFolder]);
    elseif ~statusMat
        disp(['Failed to transfer MetaData file of ' chanName ' to ' SaveFolder]);
    else
        disp(['Successful tranfer of ' chanName ' to ' SaveFolder]);
    end
end


%%%%%%%%%%%%%%%%%
outFile = [outFile, chanName];
end