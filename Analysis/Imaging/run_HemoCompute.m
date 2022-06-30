function outFile = run_HemoCompute(RawFolder,SaveFolder, varargin)
% RUN_HEMOCOMPUTE calls the function HEMOCOMPUTE from the IOI library (LabeoTech).


% Defaults:
default_Output = {'HbO.dat', 'HbR.dat'}; % This line is here just for Pipeline management.
default_opts = struct('FilterSet', 'GCaMP', 'Red', true, 'Green', true, 'Amber', true);
opts_values = struct('FilterSet',{{'GCaMP', 'jRGECO' , 'none'}},'Red',[false, true], 'Green',[false, true],'Amber',[false, true]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p, 'RawFolder', @isfolder);
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,RawFolder,SaveFolder, varargin{:});
%Initialize Variables:
RawFolder = p.Results.RawFolder;
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
outFile = default_Output;
clear p
% Translate opts to char cell array:
fields = fieldnames(opts);
fields = setdiff(fields, 'FilterSet');
idx = cellfun(@(x) opts.(x), fields);
list = fields(idx)';

% Run HemoCompute function from IOI library:
try
    HbO = HemoCompute(RawFolder,SaveFolder, opts.FilterSet, list); %#ok The output is here just to catch assertion errors from HemoCompute.
catch ME
    ME = addCause(ME, MException('umIToolbox:run_HemoCompute:UnkwnownError',...
        'Function HemoCompute Failed! Check Matlab command window for messages!'));
    throw(ME)
end
disp('Finished HemoCompute.')
end