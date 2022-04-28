function out = findMyROIfile(filename,source)
% This is a helper function for some Stats analysis functions in umIT.
% It looks for an ROI file located in the subject's folder or in a folder
% given as input.
% Inputs:
%   filename (str): Name of the ROI file to be searched.
%   source (handle): handle to an "Acquisition" or "Modality" object.
% Output:
%   out (str): full path to the ROI file.

% Append extension to filename, if not already provided:
if ~endsWith(filename, '.mat')
    filename = [filename '.mat'];
end

% Parse File path to find subject folder:
if isa(source,'Acquisition') || isa(source,'Modality')
    % If a umIT's valid object is provided, it means that the function is being run by
    % PipelineManager. In this case, find the subject's folder an try to
    % find the ROI file there:
    tmp_obj = source;
    idx = false;
    while ~idx
        tmp_obj = tmp_obj.MyParent;
        idx = isa(tmp_obj, 'Subject');
    end
    out = fullfile(tmp_obj.SaveFolder, filename);        
elseif isempty(source)
    path = fileparts(filename);
    if isempty(path)
        out = fullfile(pwd,filename);
    else
        out = filename;
    end
       else 
    error('Umitoolbox:findMyROIfile:WrongInput', 'Invalid Input. Source must be an "Acquisition" or "Modality" object');    
end

% Throw error if the file does not exist in Subject's folder:
if ~isfile(out)
    errID = 'Umitoolbox:findMyROIfile:FileNotFound';
    folder = fileparts(out); folder = strrep(folder, '\', '\\');
    errMsg = ['ROI file not found in ' folder];
    ME = MException(errID, errMsg);
    throwAsCaller(ME);
end

end
