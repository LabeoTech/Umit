function out = findMyROIfile(filename,source)
% This is a helper function for some Stats analysis functions in umIT.
% It looks for an ROI file located in the subject's folder or in a folder
% given as input.
% Inputs:
%   filename (str): Name of the ROI file to be searched.
%   source (handle): handle to an "Acquisition" or "Modality" object.
% Output:
%   out (str): full path to the ROI file.

% % Append extension to filename, if not already provided:
% [folder,filename,~] = fileparts(filename);
% filename = [filename '.mat'];
% Get protocol handle:
protObj = source;
if ~isempty(protObj)
    while 1
        protObj = protObj.MyParent;
        if ~isprop(protObj, 'MyParent')
            break
        end
    end
end
% Parse File path to find subject folder:
if isempty(protObj) 
    if isempty(folder)
        % Look for the file in the current directory
        out = fullfile(pwd,filename);
    else
        % Use fullpath from filename
        out = fullfile(folder,filename);
    end    
elseif protObj.b_isDummy
    % If the protocol object is "dummy", this means that it was created by
    % DataViewer. In this case, look for the ROI file inside the object's
    % save folder:
    out = fullfile(source.SaveFolder, filename);    
else
    
    % If the protocol object is not "dummy", look for the ROI file inside
    % the Subject's folder.
    tmp_obj = source;    
    while 1
        if isa(tmp_obj, 'Subject')
            break
        end
        tmp_obj = tmp_obj.MyParent;        
    end
    out = fullfile(tmp_obj.SaveFolder, filename);
end
clear tmp_obj protObj
% Throw error if the file does not exist in Subject's folder:
if ~isfile(out)
    errID = 'Umitoolbox:findMyROIfile:FileNotFound';
    folder = fileparts(out); folder = strrep(folder, '\', '\\');
    errMsg = ['ROI file not found in ' folder];
    ME = MException(errID, errMsg);
    throwAsCaller(ME);
end

end