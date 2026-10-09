function teardownIsolatedActiveRigFixture(fixture)
%TEARDOWNISOLATEDACTIVERIGFIXTURE Undo setupIsolatedActiveRigFixture.
%
%   teardownIsolatedActiveRigFixture(fixture)
%
%   Restores whatever Rig was Active before the fixture was created (via
%   deactivateRigTemporarily) and removes the fixture Rig's folder.

if ~isstruct(fixture)
    return
end

deactivateRigTemporarily(fixture);

if isfield(fixture, 'RigRoot') && ~isempty(fixture.RigRoot) && ...
        isfolder(fixture.RigRoot)
    rmdir(fixture.RigRoot, 's');
end
end
