function mergeRecordings(SaveFilename,folderList,filename, varargin)
% MERGERECORDINGS concatenates Image time series (Y,X,T)
% stored in a .dat file ("filename") located across the folders listed in
% "folderList".
% Inputs:
%   SaveFilename(char): full path of the ".dat" file with the merged data.
%   folderList (cell): array of full paths of the folders containing the
%       files to be merged.
%   filename (char): name of the file to be merged.
%   merge_order (int): array of integers with indices of "folderList". The data will be
%       merged in this order.


%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p, 'SaveFile', @(x) ischar(x) & ~isempty(x));
addRequired(p, 'folderList', @(x) iscell(x) & ~isempty(x));
addRequired(p, 'filename', @(x) ischar(x) & ~isempty(x));
addOptional(p, 'merge_order',[], @isnumeric);
% Parse inputs:
parse(p,SaveFilename,folderList,filename,varargin{:})
% Set merge_order variable and force it to an integer:
merge_order = p.Results.merge_order;
if isempty(merge_order)
    % If not provided, the order of merging will be the order of
    % "folderList":
    merge_order = 1:length(folderList);
else
    merge_order = round(merge_order);
end
% Append ".dat" to filenames if not done yet:
if ~endsWith(SaveFilename,'.dat')
    SaveFilename = [SaveFilename '.dat'];
end
if ~endsWith(filename,'.dat')
    filename = [filename '.dat'];
end
metaData_filename = strrep(filename, '.dat', '.mat');
% Further input validation:
% Check if merge_order contains all the indices of folderList:
assert(isequal(sort(merge_order),1:length(folderList)),'umIToolbox:mergeRecordings:MissingInput',...
    'The merge order list is incompatible with the list of folders');
% Reorder folderList following merge_order:
folderList = folderList(merge_order);
% check if the file exists in all folders:
idx = cellfun(@isfile, fullfile(folderList,filename));
if all(~idx)
    error(['The file ' filename ' was not found in any of the folders provided!']);    
elseif ~all(idx)
    disp(repmat('-',1,100))
    warning('The following folders do not contain the file %s and will be ignored:\n%s\n',...
        filename, folderList{~idx});    
    disp(repmat('-',1,100))
    % Update folderList:
    folderList = folderList(idx);    
end
% Check if the associated .mat file exists in the folder:
idx = cellfun(@(x) isfile(x), fullfile(folderList,metaData_filename));
assert(all(idx), 'umIToolbox:mergeRecordings:MissingInput',...
    'One or more associated .mat file are missing!')
% Check if the input files are image time series with dimensions {Y,X,T}:
idxDim = false(size(folderList));
idxSz = idxDim;
refSz = load(fullfile(folderList{1},metaData_filename), 'datSize');

for i = 1:length(folderList)
    md = matfile(fullfile(folderList{i}, metaData_filename));
    idxDim(i) = all(ismember(md.dim_names, {'Y','X','T'}));
    idxSz(i) = isequaln(md.datSize,refSz.datSize); 
end
assert(all(idxDim), 'umIToolbox:mergeRecordings:WrongInput',...
    'This function accepts only image time series with dimensions {"Y", "X","T"}!');
% Check if all data have the same Y,X dimensions sizes:
assert(all(idxSz), 'umIToolbox:mergeRecordings:WrongInput','All input files must have the same Y,X sizes!');
% Open new .dat file:
fidOut = fopen(SaveFilename,'w');
% Open new .mat (meta data) file:
mOut = matfile(strrep(SaveFilename, '.dat', '.mat'), 'Writable', true);
% Get some meta data from the first file:
mIn = matfile(fullfile(folderList{1},metaData_filename));
mOut.Freq = mIn.Freq;
mOut.Stim = [];
mOut.datName = mIn.datName;
mOut.datLength = [];
mOut.FirstDim = mIn.FirstDim;
mOut.dim_names = mIn.dim_names;
mOut.datFile = SaveFilename;
mOut.Datatype = mIn.Datatype;
mOut.datSize = mIn.datSize;
w = waitbar(0,'Merging data...');
if isfile(SaveFilename)
    delete(SaveFilename);
end
for i = 1:length(folderList)
    % Merge data:
    fidIn = fopen(fullfile(folderList{i}, filename),'r');
    mIn = matfile(fullfile(folderList{i}, metaData_filename));
    data = fread(fidIn,inf,['*' mIn.Datatype]);
    fwrite(fidOut,data, mIn.Datatype);
    fclose(fidIn);
    % Update meta data:
    mOut.datLength = sum([mOut.datLength, mIn.datLength]);
    mOut.Stim = [mOut.Stim, mIn.Stim];
    waitbar(i/length(folderList),w);   
end
fclose(fidOut);
close(w)
disp('Done')
end