function fixture = setupIsolatedActiveRigFixture()
%SETUPISOLATEDACTIVERIGFIXTURE Guarantee a resolvable Active/Default Rig.
%
%   fixture = setupIsolatedActiveRigFixture()
%
%   DFR-20260819-010: UMITRigStore.getActiveRig/getOrCreateDefaultRig
%   require exactly one Active Rig and error otherwise. Ambient state on a
%   given machine, or leftover corruption from an earlier test in the same
%   batch, can leave that invariant broken (zero Active Rigs). This creates
%   a throwaway Rig and makes it the sole Active/Default Rig regardless of
%   the ambient UMITRigStore state (via activateRigTemporarily), so callers
%   depending on Rig auto-resolution are hermetic. Pair with
%   teardownIsolatedActiveRigFixture(fixture) in TestMethodTeardown, which
%   restores exactly whatever this function changed and nothing more.

suffix = strrep(char(java.util.UUID.randomUUID()), '-', '');
rigStore = UMITRigStore.create(struct('rigID', ['DFR20260819010_' suffix]));

fixture = activateRigTemporarily(rigStore.getRigInfo().uuid);
fixture.RigRoot = rigStore.RigRoot;
end
