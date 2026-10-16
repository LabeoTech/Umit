classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestUmitRootIsolation < matlab.unittest.TestCase
%TESTUMITROOTISOLATION Tests never touch the real UMIT stores.
%
%   Project redesign Phase 6 (audit item 24, DFR-20261005-001):
%   getUmitFolder honours the UMIT_ROOT environment variable, and every test
%   class declares UMITTemporaryRootFixture as a shared fixture, so the
%   project store, the rig library and the preferences of the developer's
%   machine are never read or written by a test.

    methods (Test)
        function testEveryTestClassDeclaresTheFixture(testCase)
            testRoot = fileparts(mfilename('fullpath'));
            files = dir(fullfile(testRoot, '**', '*.m'));
            missing = {};
            for k = 1:numel(files)
                file = fullfile(files(k).folder, files(k).name);
                text = fileread(file);
                header = regexp(text, '^classdef[^\n]*', 'match', 'once', 'lineanchors');
                if isempty(header) || ~contains(header, 'matlab.unittest.TestCase')
                    continue
                end
                if ~contains(header, 'UMITTemporaryRootFixture')
                    missing{end + 1} = file; %#ok<AGROW>
                end
            end
            testCase.verifyEmpty(missing, sprintf(['These test classes do not ' ...
                'declare SharedTestFixtures = {UMITTemporaryRootFixture}:\n%s'], ...
                strjoin(missing, newline)));
        end

        function testStoresResolveBelowTheTemporaryRoot(testCase)
            root = getenv('UMIT_ROOT');
            testCase.assertNotEmpty(root);
            testCase.verifyEqual(getUmitFolder(), root);
            testCase.verifyTrue(startsWith(UMITProjectStore.getProjectsRoot(), root));
            testCase.verifyTrue(startsWith(UMITRigStore.getRigsRoot(), root));
            testCase.verifyFalse(startsWith(root, iDefaultRoot()));
        end

        function testEmptyOverrideMeansTheDefaultRoot(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_ROOT', ''));
            % create=false: never create the real root; '' when it is absent.
            folder = getUmitFolder('', 'create', false);
            testCase.verifyTrue(isempty(folder) || strcmp(folder, iDefaultRoot()), folder);
        end

        function testRelativeOverrideIsRefused(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_ROOT', 'relative/folder'));
            testCase.verifyError(@() getUmitFolder(), ...
                'Umitoolbox:UmitFolder:invalidRootOverride');
        end

        function testOverrideIsCreatedOnDemand(testCase)
            parent = tempname;
            testCase.addTeardown(@() iRemove(parent));
            root = fullfile(parent, 'umitRoot');
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'UMIT_ROOT', root));
            testCase.verifyEqual(getUmitFolder('rigs', 'create', false), '');
            testCase.verifyEqual(getUmitFolder('rigs'), fullfile(root, 'rigs'));
            testCase.verifyTrue(isfolder(fullfile(root, 'rigs')));
        end
    end
end

function root = iDefaultRoot()
if ispc
    home = getenv('USERPROFILE');
else
    home = getenv('HOME');
end
root = fullfile(home, 'Documents', 'LabeoTech', 'umIT');
end

function iRemove(folder)
if isfolder(folder)
    rmdir(folder, 's');
end
end
