function tmpS = protocolFcn_template(obj)
% User-defined function. This Function is used inside "Protocol" object to generate the
% list of Subjects, Acquisitions and Modalities. This template can be used
% to read recordings from LabeoTech's optical imaging systems.

% This function will be used to update the content of the obj.MainDir (Folder where the Data is).
% Input:
% obj (Protocol class) = self referring handle of Protocol class.
% Output:
% tmpS (Subject class) = array of objects from Subject class containing Acquisition and
% Modalities.

% This is a template for creating a protocol function.
% This template reads mesoscale imaging datasets created with the imaging
% systems from Labeo Technologies Inc. The imaging raw data are stored in a
% series of binary files (.bin). This template function scans a directory
% (obj.MainDir) containing folders with the imaging data and a text file
% with the recording information.
% The folder tree is organized as :
%      MainDir --
%                |
%                |--RECORDINGFOLDER-- (ex. "M2D_RS_20210401")
%                                                   |
%                                                   |-xxxx.bin
%                                                   |-xxxx.bin
%                                                   |-xxxx.bin
%                                                      ...
%                                                   |-xxxx.bin
%                                                   |-info.txt

% The names of subjects (Subject ID) and acquisitions (Acquisition ID) are encoded in the
% folder names separated by an underscore as SUBJID "_" ACQID. For example:
% Lets take the example folder above (M2D_RS_20210401), here we have:
% Subject ID: "M2D"
% Acquisition ID: "RS_20210401

% Minimal information required:
%   In order to make this function work, we need to set the following properties
%   for each object created in this function:
%       Subject:
%               -ID (str)
%       Acquisition:
%               -ID (str)
%               -Start_datetime (datetime): Timestamp of the beginning of
%               the recording.
%       Modality:
%               -ID (str)
%               -RawFolder (str): Full path for an existing folder
%                containing raw data from a given modality.
%               -RawFiles (cell): list of raw file names with the file extensions.
% * Other properties can be set as well. Please, check the documentation
% for the classes "Subject", "Acquisition", "Modality" and those derived
% from it.

% Regular expression to find the raw data files:
expLabeo = '(ai|img)_\d*.bin'; % For Labeo Data
% Look for all folders containing .bin files:
FolderNames = dir(fullfile(obj.MainDir, '**','*.bin'));
FolderNames = unique({FolderNames.folder});
% Keep only folders that start with an uppercase "M" letter:
idx = false(size(FolderNames));
recNames = cell(size(FolderNames));
for i = 1:numel(FolderNames)
    str = strsplit(FolderNames{i}, filesep);
    recNames(i) = str(end);
    if startsWith(recNames{i},'M')
        idx(i) = true;
    end
end
FolderNames = FolderNames(idx);
recNames = recNames(idx);
% Parse the list of recording folders (recNames) to obtain:
%   - Subject ID (characters between the first and second underscores)
%   - Acquisition ID (characters between the second and last underscores)
sID = cell(size(recNames)); acqID = sID;
for i = 1:numel(recNames)
    str = strsplit(recNames{i}, '_');
    sID{i} = str{1}; % Subject ID
    % rebuild Acquisition ID using all strings after subject ID.
    acqID{i} = strjoin(str(2:end), '_');
end
% Here, we create an array of subjects based on the IDs extracted from the
% folder names:
uniqS = unique(sID); % Get the list of all subjects.
% Create an array of empty Subject objects:
tmpS = Subject.empty(0,numel(uniqS));
% Loop across all folders and populate Subject objects with with Acquisition and modality objects:
for i = 1:length(uniqS)
    % Index of the current Subject:
    indx_s = find(strcmp(uniqS{i}, sID));
    % Add Subject ID and Group ID to current subject object:
    tmpS(i) = Subject(uniqS{i},'def', []);
    
    % Add Aquisitions to the current Subject:
    for j = 1:length(indx_s)
        % Create Acquisition object to the current subject:
        tmpA = Acquisition(acqID{indx_s(j)}, []);
        % Add empty FluorescenceImaging object to the Acquisition object:
        tmpA.Array.addObj(FluorescenceImaging());
        folder = FolderNames{indx_s(j)};
        % Get info for creation of FluorescenceImaging object:
        recSys = 'LabeoTech-OpticalImaging-System';
        files = getNamesFromDir(folder, expLabeo, 0, 'match');
        
        % Add Imaging modality info:
        tmpA.Array.ObjList(end).ID = 'CtxImg'; % Modality ID
        tmpA.Array.ObjList(end).RecordingSystem = recSys;  % Name of recording system
        tmpA.Array.ObjList(end).RawFiles_FP = files; % List of full paths of LabeoTech's raw imaging data.
        % Use Labeotech's "info.txt" file to gather more information about
        % the recording:
        infoFile = fullfile(folder, 'info.txt');
        if isfile(infoFile)
            txt = fileread(infoFile);
            DateTime = regexp(txt, '(?<=DateTime:\s*)\S+','match', 'once');
            DateTime = datetime(DateTime, 'InputFormat', 'yyyyMMdd_HHmmss');
            numChan = numel(regexp(txt, 'Illumination', 'match'));
            sr = str2double(regexp(txt, '(?<=FrameRateHz:\s+)\d+\.*\d*', 'match'));
            % Add info to FluorescenceImaging object:
            tmpA.Array.ObjList(end).SampleRateHz = sr; % Recording sample rate
            tmpA.Array.ObjList(end).NumberOfChannels = numChan; % number of color channels
            % Add Timestamp of the beginning of the recording to
            % Acquisition object:
            tmpA.Start_datetime = DateTime;
        end       
        % Add Acquisition Object to Subject
        tmpS(i).Array.addObj(tmpA);
    end
    
end

end