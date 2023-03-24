function outData = dummyFunc4Testing_5(data, varargin)
% FUNCTEMPLATE is a dummy function that serves as a template for creating
% functions compatible with the toolbox.
%
% The arguments have to be stated in the following order:
% 1 - Input : keywords are {File, RawFolder or SaveFolder}
% 2 - SaveIn: save directory (fullpath).
% 3 - Output: Output file name or cell array of file names.
% 4 - opts: structure containing optional parameters for the function.
%
default_Output = 'dummyFile5.dat';
default_opts = struct('ParamIntRange', 1,'ParamStr', 'val','ParamBool', false);
opts_values = struct('ParamIntRange', [1:5],'ParamStr',{{'val'}},'ParamBool',[false, true]);
% default_opts = struct('ParamIntRange', 1,'ParamNum',1,'ParamPosNum',1, 'ParamStr', 'val', 'ParamBool', false, 'ParamMultiChoice','Option1', 'ParamSingleChoice','Option1');
% opts_values = struct('ParamIntRange', [1:5],'ParamNum',[-Inf,Inf],'ParamPosNum',[eps, Inf], 'ParamStr',{{'val'}},'ParamBool',[false, true], 'ParamMultiChoice', {{'Option1', 'Option2','Option3'}'},'ParamSingleChoice',{{'Option1', 'Option2','Option3'}});%#ok  % This is here only as a reference for PIPELINEMANAGER.m.
%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
addRequired(p, 'data');
% Parse inputs:
parse(p,data, varargin{:});
%Initialize Variables:
data= p.Results.data;
%%%
% Run your code here:
disp('This is a dummy function for Pipeline testing!');
a = zeros(3,3, 'single');
save2Dat(fullfile(SaveFolder, Output), a);
out = Output;
end
