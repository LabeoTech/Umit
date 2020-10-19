classdef Protocol < handle
    % This class creates and manages a list of "Subject" objects.
    %   This class uses an user-defined function (ProtoFunc) to access the
    %   acquisition information from each subject and assign it to the
    %   pertinent objects. For more information, see documentation of the
    %   following classes: Subject, Acquisition and children of the abstract
    %   class "Modality".
    
    properties 
        MainDir % Folder containing all experiment files.
        SaveDir % Folder to save "Protocol" object and HDF5 files. Default value = current folder.
        ProtoFunc % Function handle of the user-defined OpenProtocol
        % where the Subjects and Acquisition data are created.
        Array % List of Subjects. Default: empty ObjectListManager.
    end
    properties (Access = private)
        RecordingStruct % Recording Structure obtained from ProtoFunc.
    end
    
    methods
        
        function obj = Protocol(MainDir, SaveDir, ProtoFunc, Array)
            % Class constructor.
            %   This function initiates the object "Protocol" with the
            %   properties: MainDir, SaveDir, ProtoFunc and Array. 
            %   All first inputs must be provided. If Array is empty,
            %   the function creates an emtpy Array.
            if nargin > 0
                obj.MainDir = MainDir;
                obj.SaveDir = SaveDir;
                obj.ProtoFunc = ProtoFunc;
                obj.Array = Array;
            else
                obj.Array = ObjectListManager();
            end
        end
        
        %%% Property Set Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.MainDir(obj, MainDir)
            % Set function for MainDir property.
            %   Accepts only existing Folders as input.
            mustBeFolder(MainDir); % Checks for existing Path..
            obj.MainDir = checkFolder(MainDir);
        end
        
        function set.SaveDir(obj, SaveDir)
            % Set function for SaveDir property.
            %   Accepts only existing Folders as input.
            mustBeFolder(SaveDir); % Checks for existing Path..
            obj.SaveDir = checkFolder(SaveDir);
        end
        
        function set.ProtoFunc(obj, ProtoFunc)
            % Set function for ProtoFunc property.
            %   Accepts only valid function handles.
            mustBeA(ProtoFunc, 'function_handle');
            obj.ProtoFunc = ProtoFunc; % Checks if ProtoFunc is a function handle.
        end
        
        function set.Array(obj, Array)
            % Set function for Array property.
            %   Accepts only a "ObjectListManager" object as input. If
            %   empty, creates an default "ObjectListManager" object.
            if ~isempty(Array)
                mustBeA(Array, 'ObjectListManager'); % Checks if Array is an "ObjectListManager".
                obj.Array = Array;
            else
                obj.Array = ObjectListManager();
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function generateList(obj)
           % This function uses the ProtoFunc to create the lists of
           % Subjects and Acquisitions.
         
           RecInfo = obj.ProtoFunc(obj.MainDir);
           if isempty(obj.Array.ObjList)
               genList(obj, RecInfo);
           else
               disp('There are elements in the Array. Use "Protocol".resetList to erase the Array and try again.')
           end
           
        end
        
        function resetList(obj)
            obj.Array = ObjectListManager;
        end
        
        function ignoreList = createIgnorelist(obj)
            % This function gets the full path of all files in the Protocol
            % list.
            % zu verbessern...
            ignoreList = [];
            for i = nSubj:-1:1
                nAcq = numel(obj.Array.ObjList(i).Array.ObjList);
                for j = nAcq:-1:1
                    nMod = numel(obj.Array.ObjList(i).Array.ObjList(j).Array.ObjList);
                    for k = nMod:-1:1
                        list = [];
                        folder = obj.Array.ObjList(i).Array.ObjList(j).Array.ObjList(k).Folder;
                        files = obj.Array.ObjList(i).Array.ObjList(j).Array.ObjList(k).FileName;
                        for p = numel(files):-1:1
                            list{p} = [folder files{p}];
                        end
                        ignoreList = [ignoreList; list'];                        
                    end
                end
            end
        end
        
        function delete(obj)
%             for i = 1:length(obj.Array.ObjList)
%                 obj.Array.ObjList.delete
%             end
            disp('Protocol deleted')
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% Extra Methods for Protocol class %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    end
    
end


%%% Local Functions %%%

function genList(obj, RecInfo)
% This function uses the information in the structure "RecInfo" to create
% the chain of Subject, Acquisition and Modality objects inside the "Array"
% property of the "Protocol" object.

for i = length(RecInfo):-1:1
    tmpSubj(i) = Subject(RecInfo(i).ID, RecInfo(i).GroupID, []);
    for j = length(RecInfo(i).Acquisition):-1:1
        tmpAcq(j) = Acquisition(RecInfo(i).Acquisition(j).ID,[]);
        %RecInfo(i).Acquisition.isObj = true;
        for k = 1:length(RecInfo(i).Acquisition(j).Modality)
            tmpName =  RecInfo(i).Acquisition(j).Modality(k).Type;
            eval(['tmpAcq(j).Array.addObj(' tmpName ');']);
            tmpAcq(j).Array.ObjList(k).ID = RecInfo(i).Acquisition(j).ID;
            tmpAcq(j).Array.ObjList(k).Folder =  RecInfo(i).Acquisition(j).Modality(k).Folder;
            tmpAcq(j).Array.ObjList(k).FileName = RecInfo(i).Acquisition(j).Modality(k).Files;
        end
        Mods = tmpMod;
        tmpAcq(j).Array.addObj(Mods);
    end
    Acqs = tmpAcq;
    tmpSubj(i).Array.addObj(Acqs);
end
Subjs = tmpSubj;
obj.Array.addObj(Subjs);
clear tmp*
end
   