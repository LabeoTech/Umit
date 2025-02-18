function tf = isDatStat(dataStruct)
% Validator function for .DATSTAT file.
%
% ISDATSTAT checks if the input structure has the required fields to save
% the data to a .DASTAT file.

requiredFields = {'obsID','data','b_hasEvents','b_hasMultipleMeasures','dataCategory'};
tf = all(ismember(requiredFields,fieldnames(dataStruct)));
if ~tf
    error('The input data structure is not valid!');
end

end