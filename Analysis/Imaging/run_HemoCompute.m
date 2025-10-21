function outFile = run_HemoCompute(SaveFolder, varargin)
% RUN_HEMOCOMPUTE calls the function HEMOCOMPUTE from the IOI library (LabeoTech).
% For more info, read the HEMOCOMPUTE docstring.

% Defaults:
default_Output = {'HbO.dat', 'HbR.dat'}; % This line is here just for Pipeline management.
default_opts = struct('FilterSet', genFieldsFromFile('FilterSets.mat',1), 'b_normalize', true, 'Red', true, 'Green', true, 'Amber', true, 'HbT_concentration_uM',100,'StO2perc',60);
opts_values = struct('FilterSet',{genFieldsFromFile('FilterSets.mat')'},'b_normalize',[true,false],'Red',[false, true], 'Green',[false, true],'Amber',[false, true],'HbT_concentration_uM',[eps,Inf],'StO2perc',[eps 100]);%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,SaveFolder, varargin{:});
%Initialize Variables:
SaveFolder = p.Results.SaveFolder;
opts = p.Results.opts;
outFile = default_Output;
clear p
% Translate opts to char cell array:
fields = fieldnames(opts);
fields = setdiff(fields, {'FilterSet', 'b_normalize','HbT_concentration_uM','StO2perc'});
idx = cellfun(@(x) opts.(x), fields);
list = fields(idx)';

% Run HemoCompute function from IOI library:
try
    HbO = HemoCompute(SaveFolder,SaveFolder, opts.FilterSet, list,...
        opts.b_normalize, opts.HbT_concentration_uM, opts.StO2perc); %#ok The output is here just to catch assertion errors from HemoCompute.
catch ME
    ME = addCause(ME, MException('umIToolbox:run_HemoCompute:UnkwnownError',...
        'Function HemoCompute Failed! Check Matlab command window for messages!'));
    throw(ME)
end
disp('Finished HemoCompute.')
end