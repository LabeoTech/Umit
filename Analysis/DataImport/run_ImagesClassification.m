function outFile = run_ImagesClassification(RawFolder, SaveFolder, varargin)
% RUN_IMAGESCLASSIFICATION calls the function
% IMAGESCLASSIFICATION from the IOI library (LabeoTech).

%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p, 'RawFolder', @isfolder)% For a folder as input
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure:
default_opts = struct('BinningSpatial', 1, 'BinningTemp', 1, 'b_SubROI', false, 'b_IgnoreStim', true);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Output File:
default_Output = {'fChan.dat', 'rChan.dat', 'gChan.dat', 'yChan.dat'};
addOptional(p, 'Output', default_Output, @(x) ischar(x) || isstring(x) || iscell(x));
parse(p, RawFolder, SaveFolder, varargin{:});
%Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
Output = p.Results.Output;
outFile = {};
%%%%
% Calls function from IOI library. Temporary for now.
% ImagesClassification(p.Results.RawFolder, p.Results.opts.BinningSpatial, p.Results.opts.BinningTemp, p.Results.opts.b_SubROI, p.Results.opts.b_IgnoreStim);
cd(RawFolder)
chanList = dir('*Chan.dat'); chanList = {chanList.name};
for i = 1:length(chanList)
    chanName = chanList{i};
    switch chanName
        case 'fChan.dat'
            MetaDataFileName = 'Data_Fluo.mat';
        case 'gChan.dat'
            MetaDataFileName = 'Data_green.mat';
        case 'rChan.dat'
            MetaDataFileName = 'Data_red.mat';
        case 'yChan.dat'
            MetaDataFileName = 'Data_yellow.mat';
    end
    a =  matfile(fullfile(RawFolder,MetaDataFileName), 'Writable', true);
    a.fileUUID = char(java.util.UUID.randomUUID);
    a.Datatype = 'single';
    a.datName = 'data';
    
    % TEMPORARY FIX.
    filePath = fullfile(SaveFolder,chanName);
    newmDfile = strrep(filePath, '.dat', '_info.mat');
    [~] = movefile(chanName, filePath);
    [~] = movefile(MetaDataFileName, newmDfile);
    %%%%%%%%%%%%%%%%%
    outFile = [outFile, chanName];
end
end