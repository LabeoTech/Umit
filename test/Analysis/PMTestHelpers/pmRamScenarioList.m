function scenarios = pmRamScenarioList(category)
%PMRAMSCENARIOLIST Enumerate PipelineManager RAM execution scenarios.
%
%   scenarios = pmRamScenarioList()
%   scenarios = pmRamScenarioList('fileCapable')
%   scenarios = pmRamScenarioList('ramOnly')
%
%   PipelineManager exposes exactly three RAM execution scenarios, formed by
%   the ramMode / ramSafePolicy pair (see PipelineManager's "RAM Management
%   Settings" property block):
%
%       'auto'               ramMode='auto'    - PM decides per edge whether a
%                                                DATA input is handed over as an
%                                                array or as a filename.
%       'ramsafe-strict'     ramMode='ramsafe', ramSafePolicy='strict'
%                                              - always hand over filenames;
%                                                refuse the DAG if any required
%                                                DATA input has supportsFile==false.
%       'ramsafe-bestEffort' ramMode='ramsafe', ramSafePolicy='bestEffort'
%                                              - always prefer filenames; warn
%                                                and rehydrate to RAM for DATA
%                                                inputs with supportsFile==false.
%
%   CATEGORY filters the list by what a function's DATA input can accept:
%       'all'         (default) - all three scenarios.
%       'fileCapable'          - all three scenarios. A function whose DATA
%                                input declares supportsFile==true is expected
%                                to run successfully in every scenario.
%       'ramOnly'              - all three scenarios. A function whose DATA
%                                input declares supportsFile==false is expected
%                                to run under 'auto', to be refused under
%                                'ramsafe-strict', and to warn-and-run under
%                                'ramsafe-bestEffort'.
%
%   The category argument does not change the returned list; it exists so that
%   test classes document which expectation set applies. Use
%   PMRAMSCENARIOEXPECTATION to obtain the expected outcome for a scenario.
%
%   See also buildPMForScenario, pmRamScenarioExpectation.

if nargin < 1 || isempty(category)
    category = 'all';
end

mustBeMember(lower(char(string(category))), {'all', 'filecapable', 'ramonly'});

scenarios = {'auto', 'ramsafe-strict', 'ramsafe-bestEffort'};
end
