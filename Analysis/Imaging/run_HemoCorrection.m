function outFile = run_HemoCorrection(File, SaveFolder,varargin)
% RUN_HEMOCORRECTION calls the function
% HEMOCORRECTION from the IOI library (LabeoTech).

default_Output = 'hemoCorr_fChan.dat';
%%% Arguments parsing and validation %%%
p = inputParser;
% Input File
addRequired(p,'File',@isfile)% For a file as input.
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure:
default_opts = struct('Red', true, 'Green', true, 'Amber', true);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,File, SaveFolder, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Variables:
File = p.Results.File;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
%%%%
% Calls function from IOI library. 
cd(SaveFolder)
[~,fluoMetaData]= mapDatFile(File);
if isempty(fluoMetaData)
    errID = 'MATLAB:UMIToolbox:FileNotFound';
    errMsg = ['Meta data file of fluo channel not found in ' SaveFolder];
    error(errID, errMsg);
else
    % Create temporary MAT file for compatibility with IOI function:
    copyfile(fluoMetaData.Properties.Source, 'Data_Fluo.mat');
end
% Translate opts to char cell array:
fields = fieldnames(opts);
idx = cellfun(@(x) opts.(x), fields);
list = fields(idx);
% Run HemoCorrection function from IOI library:
disp('Performing hemodynamic correction in fluo channel...')
data = HemoCorrection(SaveFolder, list);
disp('Finished hemodynamic correction.')
% delete temporary MAT file:
delete('Data_Fluo.mat');

% Save to .DAT file and create .MAT file with metaData:
datFile = fullfile(SaveFolder, default_Output);
save2Dat(datFile, data);

% Add fluo channel metaData info to new file metadata:
[~,metaData] = mapDatFile(datFile);
metaData.Properties.Writable = true;
props = setdiff(properties(fluoMetaData), properties(metaData));
for k = 1:length(props)
    eval(['metaData.' props{k} '= fluoMetaData.' props{k} ';'])
end
outFile = default_Output;
end