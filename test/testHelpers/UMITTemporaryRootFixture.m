classdef UMITTemporaryRootFixture < matlab.unittest.fixtures.Fixture
%UMITTEMPORARYROOTFIXTURE Redirect the UMIT root to a temporary folder.
%
%   Every test class declares it as a shared fixture:
%
%       classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) ...
%               TestX < matlab.unittest.TestCase
%
%   so that no test reads or writes the developer's real project store,
%   rig library or preferences (project redesign Phase 6, audit item 24;
%   DFR-20261005-001). Setup creates an empty folder and points the
%   UMIT_ROOT environment variable (honoured by getUmitFolder) at it;
%   teardown restores the previous value and deletes the folder.
%   test/TestUmitRootIsolation.m checks that every test class declares it.
%
%   Fixtures are shared between consecutive classes when their Tag is the
%   same. A class that needs a root of its own (for example an empty rig
%   library) passes a Tag:
%
%       SharedTestFixtures = {UMITTemporaryRootFixture('emptyRigLibrary')}

    properties (SetAccess = immutable)
        % Tag - Fixtures with different tags are never shared.
        Tag (1,:) char = ''
    end

    properties (SetAccess = private)
        % Root - The temporary UMIT root while the fixture is set up.
        Root (1,:) char = ''
    end

    properties (Access = private)
        PreviousValue (1,:) char = ''
    end

    properties (Constant)
        VariableName = 'UMIT_ROOT'
    end

    methods
        function fixture = UMITTemporaryRootFixture(tag)
            if nargin > 0
                fixture.Tag = char(tag);
            end
        end

        function setup(fixture)
            fixture.PreviousValue = getenv(fixture.VariableName);
            fixture.Root = tempname;
            mkdir(fixture.Root);
            setenv(fixture.VariableName, fixture.Root);
            fixture.SetupDescription = sprintf( ...
                'UMIT root redirected to the temporary folder %s.', fixture.Root);
        end

        function teardown(fixture)
            setenv(fixture.VariableName, fixture.PreviousValue);
            if ~isempty(fixture.Root) && isfolder(fixture.Root)
                try
                    rmdir(fixture.Root, 's');
                catch
                    % A file still open by a failed test must not fail the run.
                end
            end
            fixture.TeardownDescription = 'UMIT root restored.';
        end
    end

    methods (Access = protected)
        function tf = isCompatible(fixture, other)
            tf = strcmp(fixture.Tag, other.Tag);
        end
    end
end
