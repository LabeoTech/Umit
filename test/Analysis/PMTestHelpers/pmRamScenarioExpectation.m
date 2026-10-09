function expectation = pmRamScenarioExpectation(scenario, supportsFile)
%PMRAMSCENARIOEXPECTATION Expected PipelineManager outcome for a RAM scenario.
%
%   expectation = pmRamScenarioExpectation(scenario, supportsFile)
%
%   Given a scenario name from PMRAMSCENARIOLIST and whether the function's
%   DATA input declares supportsFile==true, return the outcome PipelineManager
%   is contracted to produce:
%
%       .outcome    - 'success'         : executePipeline must complete.
%                     'error'           : executePipeline must raise .identifier.
%                     'warnAndSucceed'  : executePipeline must raise .identifier
%                                         as a warning and still complete.
%       .identifier - Expected error/warning identifier, or '' for 'success'.
%
%   Contract source (PipelineManager.m):
%       validateRamSafeCompatibility - raises
%           PipelineManager:RAMSafe:UnsupportedFileInput as an error under
%           'strict' and as a warning under 'bestEffort', for any stream node
%           whose isData input has supportsFile==false. Reached from
%           executePipeline, not from validateNodePreflight.
%
%   See also pmRamScenarioList, buildPMForScenario.

arguments
    scenario (1,:) char
    supportsFile (1,1) logical
end

switch lower(scenario)
    case 'auto'
        % Auto never refuses a node on supportsFile grounds; it loads the
        % input into RAM when the port cannot take a filename.
        expectation = struct('outcome', 'success', 'identifier', '');

    case 'ramsafe-strict'
        if supportsFile
            expectation = struct('outcome', 'success', 'identifier', '');
        else
            expectation = struct('outcome', 'error', ...
                'identifier', 'PipelineManager:RAMSafe:UnsupportedFileInput');
        end

    case 'ramsafe-besteffort'
        if supportsFile
            expectation = struct('outcome', 'success', 'identifier', '');
        else
            expectation = struct('outcome', 'warnAndSucceed', ...
                'identifier', 'PipelineManager:RAMSafe:UnsupportedFileInput');
        end

    otherwise
        error('Umitoolbox:pmRamScenarioExpectation:UnknownScenario', ...
            'Unknown RAM scenario "%s".', scenario);
end
end
