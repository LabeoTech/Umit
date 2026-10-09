function text = pmFormatFileList(fileNames)
%PMFORMATFILELIST Render a file-name list for assertion messages.
%
%   text = pmFormatFileList(fileNames)
%
%   Returns a comma-separated list, or '<none>' for an empty input. Used by
%   the PipelineManager scenario tests so a failure message shows which files
%   each RAM scenario actually produced instead of only reporting that two
%   sets differed.
%
%   See also pmCollectScenarioOutputs.

if isempty(fileNames)
    text = '<none>';
    return
end

text = strjoin(cellstr(string(fileNames)), ', ');
end
