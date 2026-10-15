# Implementation notes — folder-bound rig ledger

Running log for `_tasks/completed/task_redesign_rig_management.md` (accepted 2026-10-15). Context: `AUDIT_rig_management.md`.

## Phase status

| Phase | Status |
|---|---|
| 1. Ledger core | Done, accepted (2026-10-14) |
| 2. Rig store hardening | Done, accepted (2026-10-14) |
| 3. Import paths | Done, accepted (2026-10-14) |
| 4. Analysis consumers | Done, accepted (2026-10-15) |
| 5. DataViewer | Done, accepted (2026-10-15) |
| 6. Cleanup and boundaries | Done, accepted (2026-10-15) |

All App Designer edits were applied, manually validated and promoted on 2026-10-15; the paste folder is deleted. User docs follow in `_tasks/active/task_docs_rig_ledger.md`.

## Decisions confirmed with Bruno (brief §4, 2026-10-14)

| # | Question | Decision |
|---|---|---|
| Q1 | Dual-camera import with no active rig, or no active coregistration resource | Import **uncoregistered** and write a `noRig` ledger, with a **clear warning** that asks the user to set an active coregistration. *Deviation from §1.7 "errors, never warnings"* (applies to Phase 3). |
| Q2 | Ledger that fails integrity validation | All rig-dependent steps error (`Umitoolbox:RigLedger:invalidLedger`). **No repair tool** → DFR-20261014-002. |
| Q3 | Single-camera re-import into a bound folder with committed optics | **Keep the binding**: write `reimport` only (resets the coregistration commitment). |
| Q4 | `user` / `host` in events | **Neither** is recorded, in any event (including `acknowledgeConflict`) or in the ledger lock file. *Deviation from §1.7 and §2* ("acknowledgeConflict records user and timestamp"; "all events carry user, host"). |
| Q5 | Immutable filter sets in the library | **Ledger guard only**; library versioning → DFR-20261014-003. |

## Phase 1 — Ledger core

### Files
- `RigManagement/UMITRigLedger.m` (new; static-method class, no state)
- `test/RigManagement/TestUMITRigLedger.m` (new; 31 tests)

### Public API
- Read-only: `filePath`, `exists`, `read`, `validate`, `previewRebind`
- Snapshot builders: `buildOpticsSnapshot(rigStore|RigInfo)`, `buildCoregistrationSnapshot(rigStore)` (`[]` when the rig has no active resource; the file checksum is verified against the registry), `coregistrationSnapshotFromTransform(resourceUUID, checksum, tform)`, `buildFilterSetSnapshot(id)`, `sourceRigFromRigInfo`, `bindingFromRig(rigStore, 'IncludeCoregistration', tf)`, `makeBinding(sourceRig, optics, coreg)`
- Hashing: `canonicalJSON`, `computeHash`, `shortHash`
- Writers (all require `'Origin'`): `writeNoRig`, `bind`, `rebind`, `revertBinding`, `commit`, `recordReimport`, `acknowledgeConflict`, `clearStaleLock`
- Guards: `assertCanUseOptics`, `assertCanUseCoregistration('RequestedHash', h)`, `assertCanUseFilterSet(folder, snapshot)`
- Shared error text: `buildErrorMessage(operation, folder, why, fixes)`

### Deviations from the brief and design choices
1. **Snapshots are stored as lists, not hash-keyed objects.** §2 shows `"optics": { "<hash>": {...} }`. MATLAB `jsondecode` turns keys into valid field names (`sha256:` → `sha256_`) and truncates them at 63 characters (`namelengthmax`), so a 71-character key cannot survive a round trip. Each snapshot is an object in a list and carries its own `hash` field, which is the key. On read, every snapshot's recomputed hash must equal its `hash`, and duplicates are invalid.
2. **No `user` / `host`** (Q4). The base event fields are `eventUUID, type, timestamp, toolboxVersion, origin`.
3. **`toolboxVersion`** is the constant `UMITRigLedger.ToolboxVersion = '2.0.0'`, mirroring `<param.version>` in `Umit.prj`. The repository has no runtime version API, and reading the `.prj` file would break in the compiled DataViewer. *Open issue:* keep it in sync manually, or add a shared version function later.
4. **What acknowledgement means.** `acknowledgeConflict` accepts the current binding: each conflicting part is re-committed to the currently bound hash (replay sets `commitments.<part> = binding.<part>Hash`). Without this, the folder would stay Conflicted forever. Outputs computed before the rebind are not recomputed, and `hadConflict` stays true. The event records `parts`, `acceptedOpticsHash`, `acceptedCoregistrationHash` and `comment`.
5. **Revert (conflict fix 1)** is `revertBinding`: a rebind to the committed snapshots stored in the ledger. It needs no library rig. Its `sourceRig` is the rig that originally introduced the committed optics hash.
6. **How a conflict ends.** "Conflicted" is derived: the folder is bound and a committed optics/coregistration hash differs from the bound hash. It ends by a revert, by an acknowledgement, or, for a coregistration-only conflict, by a `reimport` (which resets the coregistration commitment).
7. **Rebind to the same hashes and the same source-rig UUID writes nothing.** `bind` on a bound folder throws `alreadyBound` (use `rebind`). `rebind` on an Unknown/NoRig folder throws `notBound`.
8. **`writeNoRig`**: Unknown → creates the ledger; NoRig → no-op (returns false); bound → `alreadyBound`. Once a folder is bound it can never go back to NoRig (Q3 keeps bindings).
9. **Coregistration guard.** `RequestedHash` must equal the bound coregistration. A different hash with a committed coregistration throws `coregistrationConflict` (a different transform can only come through re-import). A different hash without a commitment throws `coregistrationNotBound` (rebind first; this is free).
10. **Filter-set commits require a bound, non-conflicted folder** (NoRig never commits). Filter sets are not part of the binding, so they can never cause a rebind conflict. Only rule §1.5 applies.
11. **Records without a spectrumID** (camera/illumination) are kept in the optics snapshot with `spectrumID ''`, `spectrumSource ''` and an empty `response`. Consumers that need them will error in Phase 4.
12. **`spectrumSource`** is `'builtIn'` when the spectrum file comes from `RigManagement/resources`, `'user'` otherwise.
13. **Locking.** A folder lock `UMITRigLedger.lock/` in the save folder, created atomically with `java.io.File.mkdir` (same approach as `UMITRigStore.iAcquireWriteLock`). It contains `lockInfo.mat` (operation, createdOn, processID; no user/host). A held lock throws `ledgerLocked` and reports the lock age. `clearStaleLock(folder, 'MinimumAge', minutes(10))` removes a leftover lock.
14. **Atomic write.** The writer encodes to a temp file in the save folder, re-reads it and runs the full integrity check plus replay on the candidate, moves the old file to a backup, renames the candidate into place (restoring the backup on failure), then deletes the backup. This is the `saveMatAtomic` pattern. The whole file is rewritten on each event; the event list is append-only.
15. **Reads never write.** `read`, `validate`, `previewRebind` and the guards only read. A folder without a ledger reads as `Unknown` and no file is created (tested).

### Canonical hash form (normative; also in the class help)
- `hash = 'sha256:' + lowercase hex(SHA-256(UTF-8 bytes of canonicalJSON(snapshot)))`, via Java `MessageDigest`.
- The top-level `hash` field is excluded.
- Objects: keys sorted by ordinal (code-point) order; `{"key":value,...}`; no whitespace anywhere.
- Lists: struct arrays, cell arrays and numeric vectors become `[v1,v2,...]` whatever their orientation. A numeric matrix becomes a list of row lists (e.g. the 3x3 transform). The fields `cameras` and `illuminations` are always lists, even with one element.
- Numbers: `0` → `0`. Any other finite value is rounded to 6 significant digits and written as `<mantissa>e<exponent>`, with trailing mantissa zeros removed, no `.` for integer mantissas, and a plain signed integer exponent. Examples: `1e0`, `5e-1`, `-1.23457e6`. NaN/Inf are rejected.
- Text: `"..."` with `\"`, `\\` and control characters escaped as `\u00XX`. Logical scalars: `true`/`false`. Empty numeric: `[]`.
- Builders round every stored number to 6 significant digits before storage. Re-rounding is then idempotent and the hash survives the `jsonencode`/`jsondecode` round trip. (So a ledger snapshot's spectra are 6-significant-digit copies of the library spectra.)
- Snapshot content: hashes ignore rig identity (two rigs with identical optics share an optics hash). Cameras are sorted by `index` and illuminations by `name`, so the library's record order does not matter.

### Snapshot schemas
- **optics**: `wavelengthNm[301]`, `cameras[{index, spectrumID, model, spectrumSource, response[301]}]`, `illuminations[{name, spectrumID, spectrumSource, response[301]}]`, `hash`
- **coregistration**: `resourceUUID`, `sourceChecksum` (registry SHA-256), `transformClass` (e.g. `affine2d`), `transformMatrix` (3x3; `T` for `affine2d`-family objects, `A` for `*tform2d` objects, as stored by its class), `hash`
- **filterSet**: `id`, `excitationSpectrumID`, `emissionSpectrumID`, `excitationResponse[301]`, `emissionResponse[301]` (an empty or `none` spectrum is all ones), `hash`

### Event schemas
Base: `eventUUID, type, timestamp (ISO-8601 UTC), toolboxVersion, origin`.
- `noRig`
- `bind`: `sourceRig{uuid, rigID, displayName, libraryModifiedOn}`, `opticsHash`, `coregistrationHash` (null when there is none)
- `rebind`: as `bind`, plus `previousOpticsHash`, `previousCoregistrationHash`
- `commit`: `part`, `hash`, `filterSetID` (null unless filterSet), `firstUsedBy` (defaults to `origin`)
- `reimport`: `previousCoregistrationHash` (informational)
- `conflict`: `parts`, `committedOpticsHash`, `committedCoregistrationHash`, `requestedOpticsHash`, `requestedCoregistrationHash`; it must directly follow the rebind that caused it
- `acknowledgeConflict`: `comment`, `parts`, `acceptedOpticsHash`, `acceptedCoregistrationHash`

### Integrity rules (`validate`)
- `format`, `schemaVersion == 1`, `ledgerUUID`, `snapshots.{optics,coregistration,filterSet}`, and a non-empty `events` list.
- Every snapshot has exactly its schema fields and its recomputed hash equals its `hash`; there are no duplicate snapshots.
- Every event has the base fields and a unique `eventUUID`.
- Event replay:
  - `noRig` only before any binding.
  - `bind` only when unbound; `rebind` only when bound, with `previous*` matching the binding it replaces.
  - Referenced snapshots exist.
  - `commit` only while bound and not Conflicted; the hash equals the bound hash; no duplicate commitment of a part or a filter-set ID.
  - `conflict` exactly after a rebind that changed a committed part.
  - `acknowledgeConflict` only while Conflicted.
  - Unknown event types are invalid.

### Error identifiers
All guard and writer refusals use `buildErrorMessage`: *"<Operation> was refused for folder "<folder>". Why: … How to fix: (1) … (2) …"*. Rig display names and short hashes (`sha256:` + 12 hex) appear wherever a committed or requested state is named. Conflict messages always offer (1) `revertBinding`, (2) a fresh save folder / re-import, and (3) `acknowledgeConflict`.

| Identifier | When |
|---|---|
| `Umitoolbox:RigLedger:noLedger` | Guard/re-import on a folder without a ledger |
| `Umitoolbox:RigLedger:noRigBound` | Guard/commit on a NoRig folder |
| `Umitoolbox:RigLedger:conflicted` | Guard/commit on a Conflicted folder |
| `Umitoolbox:RigLedger:invalidLedger` | Any read/guard/writer on a ledger failing integrity validation |
| `Umitoolbox:RigLedger:alreadyBound` | `bind` or `writeNoRig` on a bound folder |
| `Umitoolbox:RigLedger:notBound` | `rebind`/`revertBinding` on an unbound folder |
| `Umitoolbox:RigLedger:noConflict` | `acknowledgeConflict` without an open conflict |
| `Umitoolbox:RigLedger:commitMismatch` | Committing an optics/coregistration hash that is not the bound one |
| `Umitoolbox:RigLedger:noCoregistrationBound` | Coregistration guard/commit when the binding has no transform |
| `Umitoolbox:RigLedger:coregistrationConflict` | A requested transform differs from the committed one |
| `Umitoolbox:RigLedger:coregistrationNotBound` | A requested transform differs from the bound, uncommitted one |
| `Umitoolbox:RigLedger:filterSetModified` | Filter-set ID committed with different content (§1.5) |
| `Umitoolbox:RigLedger:hashMismatch` | A supplied snapshot's `hash` does not match its content |
| `Umitoolbox:RigLedger:invalidSnapshot` | A supplied snapshot has missing or unexpected fields |
| `Umitoolbox:RigLedger:invalidCoregistrationResource` | Active coregistration file fails its checksum or has no `tform` |
| `Umitoolbox:RigLedger:ledgerLocked` / `lockNotStale` / `lockCreationFailed` | Locking |
| `Umitoolbox:RigLedger:readFailed` / `writeFailed` / `hashFailed` / `invalidInput` | I/O and argument errors |

### Validation
`validateTaskChanges(ChangedFiles=[UMITRigLedger.m, TestUMITRigLedger.m], TestTargets=TestUMITRigLedger.m, Mode="targeted")` → **passed**: 31/31 tests, 0 static messages.

### Open issues (Phase 1)
- `toolboxVersion` constant (see item 3).
- The atomic-write tests cover the lock, rejected input, absence of leftover files and re-read verification. A failure of the final rename itself is not simulated, because it cannot be triggered deterministically on Windows without mocking.

## Phase 2 — Rig store hardening

### Files
- `RigManagement/UMITRigStore.m` (modified)
- `Analysis/DataImport/ImagesClassification.m` (dead `'createdDefault'` case removed from `localDescribeResolution`)
- `test/RigManagement/TestUMITRigStore.m` (modified; 6 new tests, 2 tests rewritten, 1 test adjusted)

### Changes
1. **`listRigs` guards the `defaultRig.mat` read.** A malformed pointer now reads as "no default" (`IsDefault` all false), so `validateLifecycle`, `getActiveRig`, `listAssignableRigs` and the Rig Manager list keep working. `getDefaultRig` (and `getOrCreateDefaultRig` through it) still throws `Umitoolbox:UMITRigStore:invalidDefaultRig`, so the problem stays visible.
2. **Stale store-lock recovery.**
   - New `UMITRigStore.getStoreLockInfo()` (read-only: `exists`, `operation`, `createdOn`, `age`).
   - New `UMITRigStore.clearStaleStoreLock('MinimumAge', minutes(10))`. It returns false when no lock exists and throws `Umitoolbox:UMITRigStore:storeLockNotStale` for a younger lock.
   - `iAcquireStoreLock` now creates the lock atomically (`java.io.File.mkdir`; the old `isfolder` check plus `mkdir` could race) and records `processID`.
   - The `storeLocked` message now names the blocking operation and the lock age, and points to `clearStaleStoreLock`.
3. **`iUpgradeRigInfo` no longer defaults a missing `status` to `'active'`.** A pre-v2 record without `status` loads as `'available'`, unless it carries the legacy default marker `isDefault = true`, in which case it loads as `'active'`. This keeps the legacy-default promotion in `getOrCreateDefaultRig` consistent: `iFindLegacyDefaultUUID` already treats such a record as active, and `getDefaultRig` requires `active`. Several legacy `isDefault` rigs were already an explicit error (`multipleLegacyDefaults`).
4. **No more `OiS200` auto-creation.** In `getOrCreateDefaultRig`, an empty or all-archived library now throws `Umitoolbox:UMITRigStore:noRigInLibrary`, with a message pointing to the Rig Manager template and `importRig`. Its `wasCreated` output is always false and is kept only for the callers that Phases 3–4 replace. The unused `iMakeDefaultRigID` was removed.
   - `OiS200.json` stays as a template: `getBuiltInRigDefinition('OiS200')` is unchanged, and `RigManagerTool.mlapp` already offers it (verified by grepping the app's `document.xml`; no app change needed).
   - Callers: `ImagesClassification`, `importFromTif` and `DataViewer_Coreg2Cams.mlapp` still call `getOrCreateDefaultRig`; they are rewritten in Phases 3–4. **Interim regression:** on a machine with an empty rig library, these three now error instead of silently creating `OiS200`.
5. **Rig export/import** (new):
   - `filePath = store.exportRig(filePath)` writes a self-contained `.umitrig` file. It is MAT format, with variable `UMITRigExport` holding `format = 'umit-rig-export'`, `schemaVersion = 1`, `exportedOn`, `rig` (the full RigInfo), `spectra` (category, id, metadata and 301-sample response of every camera/illumination spectrum the rig references) and `resources` (registry record plus file bytes of every managed coregistration resource, active, available or archived). The rig must pass full validation. The file is written to a temp file, re-loaded and then renamed. Filter sets are not part of a rig and are not exported.
   - `store = UMITRigStore.importRig(filePath, 'RigID', id, 'NewIdentity', tf)`:
     - Keeps the rig UUID, timestamps, resource UUIDs, checksums, statuses and the active coregistration pointer. As a result, ledger optics and coregistration snapshot hashes are identical after a round trip (tested).
     - A UUID already in the library → `importRigExists` (use `NewIdentity`); an existing rig ID → `importRigIDExists` (use `RigID`).
     - A spectrum missing from the library is written to the user library; a same-ID spectrum with different content → `importSpectrumConflict` (spectrum IDs are immutable).
     - All spectra are checked before anything is written. The rig is staged in a temp folder, resource checksums are verified, and the rig is installed under the store lock and fully validated (removed again if invalid).
     - The rig becomes Active (and the default pointer is written) only when no Active rig exists, otherwise Available. A foreign or malformed file → `invalidRigExport`; a failed export → `exportRigFailed`.
     - Caveat: spectra are written to the user library just before the rig folder is installed. If installing the rig folder then fails, those spectra stay, which is harmless because their content matches their ID.

### New error identifiers (UMITRigStore)
`noRigInLibrary`, `storeLockNotStale`, `exportRigFailed`, `invalidRigExport`, `importRigFailed`, `importRigExists`, `importRigIDExists`, `importSpectrumConflict` (all `Umitoolbox:UMITRigStore:*`).

### Tests
New in `TestUMITRigStore`:
- `testLegacyRigWithoutStatusLoadsAvailable`
- `testListRigsToleratesMalformedDefaultPointer` (backs up and restores the real `defaultRig.mat`)
- `testStaleStoreLockCanBeCleared` (skipped if a real store lock exists)
- `testExportImportRoundTripIsSelfContained`
- `testImportRestoresMissingUserSpectrumAndRejectsConflicts`
- `testExportRejectsInvalidRigAndImportRejectsForeignFile`

Changed:
- `testGetOrCreateDefaultRigCreatesWhenNoneExist` → `testGetOrCreateDefaultRigNeverCreatesARig` (expects `noRigInLibrary`).
- `testGetOrCreateDefaultRigReusesCreatedDefault` removed: there is no creation any more.
- `testV1UpgradePreservesCoregistration` now sets the legacy `isDefault` marker so the record still loads Active.
- `testGetOrCreateDefaultRigErrorsWithMultipleNonDefault` now sets up two Active rigs (see "Empty-library run" below).

### Validation
- `validateTaskChanges(ChangedFiles=[UMITRigStore.m, UMITRigLedger.m, ImagesClassification.m, TestUMITRigStore.m], TestTargets="test/RigManagement", Mode="subsystem")`: 112 run, 108 passed, **0 failed**, 4 incomplete.
- The 4 incomplete tests are the existing guarded `getOrCreateDefaultRig` tests, which only run on a machine with an empty rig library. The validator reports status "failed" because of them.
- 3 static advisories in `ImagesClassification.m` (FNDSB/LOGSUM/ISCL) were there before this task and are untouched.
- Re-run on `TestUMITRigStore` after the pragma cleanup: 53 run, 49 passed, 0 failed, 4 incomplete (same skips), 0 static messages.
- The real rig store was checked afterwards: same 16 rigs, `OiS200` still Active and default, no store lock left.

### Empty-library run (2026-10-14, at Bruno's suggestion)
The 4 guarded tests were run against a truly empty library by a scratchpad helper (`runEmptyRigLibraryTests.m`) that:
1. zips the real rigs root and checks that every file is in the zip;
2. moves the real root aside by renaming it (not deleting it);
3. runs the tests through `validateTaskChanges` with `TestMethods`;
4. always moves the real root back;
5. compares the restored file list and the `listRigs` UUID/status/default table with the originals.

Results:
- First run: 3/4 passed. `testGetOrCreateDefaultRigErrorsWithMultipleNonDefault` was **stale, not a Phase 2 regression**. It created two rigs, but `create` makes the first rig of an empty library Active, so `getOrCreateDefaultRig` correctly promoted it (a code path Phase 2 did not change). The test was rewritten to reproduce the real ambiguous state (two Active rigs, with the second record written directly).
- Second run: **4/4 passed**.
- Both runs restored the library identically (24 files, same rigs, statuses and default) with no test leftovers.
- Backup zips are kept in the session scratchpad. The same procedure will cover Phase 3's "import on a machine with an empty library" tests.

### Open question (Phase 2)
- Export/import is backend-only. Adding Export/Import buttons to `RigManagerTool.mlapp` would be a strict-mode `PASTE_INTO_MLAPP.txt` change; it is not requested by the brief and has not been done.

## Phase 3 — Import paths

### Decisions confirmed with Bruno (2026-10-14, before Phase 3)
| # | Case | Decision |
|---|---|---|
| P3-A | Dual-camera import with an active rig and transform, but `applyTform2Cams` fails | Handled **like Q1**: warn, keep the data uncoregistered, record "no rig" in a new folder; nothing is committed. |
| P3-B | Dual-camera **re-import** into a bound folder with no Active rig or no active transform | **Keep the optics binding and drop the coregistration**: write `reimport`, then rebind to the folder's own optics snapshot with no coregistration (`UMITRigLedger.dropCoregistrationBinding`), and warn. |
| P3-C | `ImagesClassification(..., 'ApplyRig', false)` (calibration source data) | **No ledger reads or writes.** |

### Files
- `Analysis/DataImport/ImagesClassification.m`: rig block rewritten; `AcqInfos.mat` saved with `saveMatAtomic` (O-12); header help updated.
- `Analysis/DataImport/importFromTif.m`: rig resolution removed; `noRig` ledger; header help updated.
- `RigManagement/UMITRigLedger.m`:
  - new `assertCanReimport` (read-only preflight, error `Umitoolbox:RigLedger:reimportConflict`);
  - new `dropCoregistrationBinding` (writer);
  - **fix**: `read` now reports an unreadable or non-JSON ledger as `invalidLedger`. Before, it surfaced the internal `readFailed`; the writers and `validate` already did this correctly. Found by the new import tests.
- Tests:
  - `TestUMITRigLedger`: 2 new tests, 1 assertion added.
  - `TestImportFromTif`: rig fixture removed, 1 test updated, 3 new tests.
  - `TestRunImagesClassification`: rig fixture removed, 2 tests updated, 1 obsolete test replaced, 4 new tests.
  - `TestRunImagesClassificationCoregistration`: 2 tests extended, 4 new tests, new helper `createCoregisteredRig`.
  - New class `TestImportEmptyRigLibrary`.

### ImagesClassification behaviour
1. Rig information goes only to the ledger. All `AcqInfos.mat` rig-field writes were removed (the file is rebuilt from `info.txt`, so legacy fields disappear on re-import).
2. Dual-camera detection uses `ReadInfoFile(RawFolder).MultiCam` (the same `info.txt` flag that classification uses), so it is known **before** classification.
3. Preflight, before `SaveFolder` is touched:
   - `UMITRigLedger.read` (an invalid ledger is an error, Q2);
   - for dual-camera data: resolve the **Active** rig (`validateLifecycle`/`getActiveRig`; several Active rigs is an error; no default pointer and no `getOrCreateDefaultRig`);
   - load and validate its active transform (an invalid resource is still `InvalidActiveCameraCoregistration`);
   - build the binding (`bindingFromRig`; checksum verified);
   - run `assertCanReimport`.
4. On a re-import (the ledger exists), `reimport` is written **before** classification overwrites files, so the coregistration commitment is never left claiming data that is being replaced.
5. After classification:
   - **Single camera**: first import → `writeNoRig`. Re-import → nothing more (Q3).
   - **Dual camera, transform applied**:
     - Unknown/NoRig folder → `bind`; bound folder → `rebind`. No conflict is possible after the preflight and the reset; one would be an `ImagesClassification:RigConflict` error.
     - Then `commit` coregistration.
     - Console report of the rig and resource used; `DataParams.cameraCoregistration` provenance is unchanged (best-effort).
   - **Dual camera, not coregistered** (no Active rig, no active transform, or application failed):
     - warning `Umitoolbox:ImagesClassification:ImportedWithoutCoregistration` (states the reason and how to fix it);
     - Unknown → `noRig`; NoRig → unchanged; bound → `dropCoregistrationBinding` (P3-B).
6. Output `rigResolution` is now a ledger summary:
   - `ledgerState`, `rigUUID`, `rigID`, `rigDisplayName`, `opticsHash`, `coregistrationHash`;
   - `cameraCoregistration.{status, wasApplied, skipReason, resourceUUID, displayName, fileName, managedFilePath, dataParamsRecorded}`;
   - `skipReason` ∈ `singleCameraAcquisition`, `noActiveRig`, `noActiveResource`, `applicationFailed`;
   - `wasCreated` and `resolution` were removed. No other code consumes this output: grepped `.m` files and the DataViewer/RigManagerTool/Coreg2Cams `document.xml`.
7. Removed: the warning `DefaultRigAutoResolved` and the `localDescribeResolution`/`localEmptyCoregistrationResolution`/`localResourceResolution` helpers. `genBackupFolder` does not manage `*.json`, so the ledger survives `ERASE`/`GENBACKUP` (verified in its file list).

### importFromTif behaviour
- No rig resolution. The ledger is read before any file is written (an invalid ledger is an error).
- After a successful import: `writeNoRig` only when the folder has no ledger. A bound or NoRig ledger is left unchanged.
- TIFF imports append channels and never overwrite raw-derived data, so they write no `reimport`.
- `AcqInfos.mat` no longer gets `rigUUID`/`rigID`, and **legacy fields are dropped** when the file is rewritten.
- `rigResolution` = `struct('ledgerState', ...)`.
- The archived-rig error (`importFromTif:ArchivedRig`) is gone with the rig resolution.

### Tests and validation
- Subsystem run: `validateTaskChanges(TestTargets=[TestUMITRigLedger, test/Analysis/importFromTif, test/Analysis/ImagesClassification, test/Analysis/DataImport], Mode="subsystem")` → 100 run, 95 passed, 2 failed, 3 incomplete.
  - The 2 failures were the `read` gap above.
  - The 3 incomplete were the empty-library guards.
  - The dual-camera tests ran on `D:\UMIT-DEV\Test2Cam` and passed.
- After the fix: `TestUMITRigLedger`, `TestImportFromTif` and `TestRunImagesClassification` → **75/75 passed**, 0 static messages.
- `TestGetEvents` (pipeline with ImagesClassification): 3/3 passed.
- Empty-library run (real library moved aside and restored, same procedure as Phase 2): `TestImportEmptyRigLibrary` **3/3 passed**. It covers single-camera → NoRig; dual-camera → uncoregistered + NoRig + warning with `noActiveRig`; TIFF → NoRig; the library stays empty. The library was restored identically (24 files, same rigs).
- Static: 3 FNDSB/LOGSUM/ISCL advisories in `ImagesClassification.m` were there before this task (untouched code).

Brief Phase 3 test list:

| Required test | Covered by |
|---|---|
| Single-camera import | `testReturnsRelativeFileManifest`, `testSingleCameraImportLeavesRigLibraryUntouched` |
| Dual-camera import | `testPrefersActiveRigResourceOverLegacyFallback`, `testFallsBackGracefullyWithNoActiveRigResource` |
| Re-import, same rig | `testReimportWithSameRigRecommitsCoregistration`, `testSingleCameraReimportKeepsBinding`, `testRepeatedSingleCameraImportRecordsReimport` |
| Re-import, different rig | `testReimportWithDifferentRigBeforeOpticsCommitRebinds`, `testReimportWithDifferentOpticsAfterCommitIsRefused` (also checks no file was touched), `testReimportWithoutActiveResourceKeepsOpticsBinding` |
| Empty library | `TestImportEmptyRigLibrary` |

### Interim behaviour until Phases 4–5 (known, intended to be fixed there)
- `DatImageSource` still calls `ensureDatasetRigAssociation` when DataViewer opens a `.dat`. New imports no longer carry `AcqInfos` rig fields, so **opening a freshly imported folder in DataViewer stamps the Active rig into its `AcqInfos.mat`** (O-1). The ledger is not affected. Fixed in Phase 5.
- `HemoCompute` still resolves optics from `AcqInfos.mat` rig fields, falling back to the Active rig with a warning. Fixed in Phase 4.
- `DataViewer_Coreg2Cams` still uses `getOrCreateDefaultRig` (errors on an empty library since Phase 2). Phase 4.

## Phase 4 — Analysis consumers

### Decisions confirmed with Bruno (2026-10-14, during Phase 4)
| # | Case | Decision |
|---|---|---|
| P4-A | Multi-folder pipeline where one folder's rig-guard preflight fails | **Skip that folder only.** It is logged as `Failed` with the guard message and the other folders run. If every folder fails, the run still errors before any work starts. |
| P4-C | Step skipping after an acknowledged conflict | **Recompute**: outputs of rig-dependent steps older than the latest `acknowledgeConflict` are never reused (see "Step skipping" below). |
| P4-D | Filter set missing from the library, committed copy used | **Dedicated warning** `Umitoolbox:HemoCompute:FilterSetFromLedger`. |
| P4-B | Which library rig `DataViewer_Coreg2Cams` saves a calibration into | **The folder's bound rig** (an error if it is archived or not in this library). Otherwise the rig passed by the caller. Otherwise a picker dialog with the Active rig pre-selected. Never the machine default. |

### Files
- `Analysis/Imaging/HemoCompute.m`: optics come from the ledger only; guards; commits; help updated.
- `RigManagement/UMITRigStore.m`:
  - new static `opticalConfigurationFromSnapshots` (pure computation over ledger snapshots; same `Umitoolbox:UMITRigStore:*` resolution error IDs as before);
  - `resolveOpticalConfiguration` now snapshots the library rig and delegates to it, so there is one implementation and `TestUMITRigOptics` keeps covering it.
- `RigManagement/UMITRigLedger.m`: new read-only `committedFilterSet(saveFolder, id)`.
- `classes/PipelineManager.m`:
  - `validateNodePreflight` adds `report.folderErrors`;
  - `executePipeline` skips and logs folders with a folder-specific preflight error;
  - new private `demoteFolderPreflightErrors`, used by `diagnosePipeline`/`autoValidate` so the GUI's pre-run check does not block a batch for one bad folder.
- New `GUI/buildLibraryRigChoices.m` (pure), `GUI/promptSelectLibraryRig.m` (modal picker), `GUI/DataViewer/resolveCoregistrationRig.m` (Coreg2Cams target resolution).
- `PASTE_INTO_MLAPP.txt`: 4 edits to `DataViewer_Coreg2Cams.mlapp` (strict mode). This replaces an older batch for the DataViewer import menu, which was already applied: the app already adds `'ImagesClassification'`.
- Tests:
  - `TestHemoCompute`: fixture binds through the ledger from a RigInfo struct (no library rig); 2 bypass tests replaced; 7 new ledger tests.
  - `TestRunHemoCompute`: fixture likewise; the Active-rig-fallback test replaced; 1 new per-folder pipeline test.
  - New `test/DataViewer/TestCoregistrationRigResolution.m`.

### HemoCompute
- Order:
  1. `UMITRigLedger.assertCanUseOptics`: no ledger, NoRig, Conflicted or invalid → error, with the shared message that explains how to bind, resolve or repair.
  2. Filter-set snapshot: from the library; if the library no longer has the ID (or one of its spectra), the folder's committed snapshot is used (`committedFilterSet`), with the dedicated warning **`Umitoolbox:HemoCompute:FilterSetFromLedger`**. This was requested by Bruno (2026-10-14). It is not a §1.7 violation: nothing is refused and the folder's committed content is used. The warning names the filter set, why the library load failed, the folder and the short hash, and says how to restore the library copy. It is raised by both the run and the preflight. Tested in `testCommittedFilterSetMissingFromLibraryIsUsed`, which also checks that a library filter set is warning-free.
  3. `assertCanUseFilterSet` (§1.5).
  4. AcqInfos channel mapping (unchanged; acquisition metadata, not rig data).
  5. `opticalConfigurationFromSnapshots`.
- `preflight` runs exactly the same resolution, so pipelines fail during preflight.
- After a successful computation and output check: commit optics and the filter set, each only when not already committed (`firstUsedBy = 'HemoCompute'`).
- Removed:
  - the Active/default-rig fallback (`localResolveActiveRig`, warning `MissingRigUsingActiveRig`, error `MissingRig`);
  - reading `AcqInfos` `rigUUID`/`rigID`;
  - **the internal `'OpticalInfo'` Name-Value**. It let callers inject optics around the ledger guards, against §1.1 ("analysis code reads rig information ONLY from the ledger"); only two tests used it. The internal `OpticalInfoChannelMismatch` assertion stays as defensive code.
- Numerics: snapshot spectra are 6-significant-digit copies. The legacy numerical regression in `TestUMITRigOptics` (AbsTol 1e-3), which now runs through the delegating store method, still passes.

### DataViewer_Coreg2Cams (strict mode, pending manual application)
- Finding: Coreg2Cams **never applies a transform to a data folder**. It classifies raw calibration data into a temporary folder (`ApplyRig=false`), computes a transform and saves it into a library rig; the `DataParams` record goes into that temporary folder. "Applying a different transform to a committed folder" therefore only happens through `ImagesClassification`, which the ledger guards (`assertCanReimport`, the reimport reset, and the coregistration commit). No conflict check is needed in Coreg2Cams; saving into the library never touches folder data (§1.1).
- Change (P4-B): `ensureRigContext` uses `resolveCoregistrationRig(app.DataFolder, 'RigUUID', app.RigUUID)`, then `promptSelectLibraryRig` when needed. `startupFcn` no longer treats a passed RigUUID as already resolved, so the folder ledger takes priority over the project-session mirror DataViewer passes (C-6). Alerts now say why no rig is available. After a save, the status tells the user that earlier imports keep their transform.
- **Pending**: Bruno applies `PASTE_INTO_MLAPP.txt`. Its 4 OLD blocks were checked to match the current app source exactly once each. Then `exportMatlabAppToDelta('DataViewer_Coreg2Cams')`, the manual checks listed in the paste file, then `confirmAndUpdateReference`.
- DataViewer itself still passes the project-session RigUUID/RigID; the ledger now takes priority over it. Removing that argument is Phase 5/6 (C-6).

### PipelineManager
- Guard errors at **run time** (inside the step) already stop only their folder: the per-folder `try/catch` logs `Failed` with the message and continues (unchanged).
- Guard errors at **preflight** (P4-A): `executePipeline` skips and logs only the failing folders. A non-folder preflight error, or all folders failing, still throws `PipelineManager:validateNodePreflight:Failed`. `diagnosePipeline`/`autoValidate` report skippable folder errors as warnings ("This folder will be skipped and logged as Failed: ...").
- No interactive prompts: the guards never prompt, HemoCompute has none, and `ImagesClassification` defaults `backupOpts` to `GENBACKUP` under PipelineManager (unchanged). The only prompts are in the GUI tool (PipelineManagerTool), not in batch execution.

### Step skipping: does it need rig-awareness? (brief §3, Phase 4)
**Conclusion: no, with one intended exception.** `isHistoryMatch` compares step name, version, parameters and input sources, not rig state. That is sufficient because:
1. Preflight runs for every folder before skip planning, whether or not the step will be skipped. A Conflicted, unbound or NoRig folder, or a modified filter set, is refused before any output can be reused.
2. In a usable folder, every committed part equals the bound part and commitments are permanent. A HemoCompute output exists only if HemoCompute succeeded, and success always commits the optics and filter set (the step only counts as executed, and gets a `dataHistory` entry, when the commit succeeded). So a reusable output was computed with exactly the optics the folder is still committed to. Filter sets: the ID is a recorded parameter, and same ID with different content is refused by the guard.
3. A re-import erases `dataHistory.mat`: `genBackupFolder` manages `*.mat`, verified. So outputs from before a re-import are never matched afterwards.

**Exception: an acknowledged conflict — now handled (decision P4-C, Bruno 2026-10-14).** `acknowledgeConflict` re-commits the folder to the new optics, but outputs computed before the rebind stay on disk with their `dataHistory` entries. With `b_skipSteps = true` (the default), a re-run would have reused them, giving a folder that says rig B with Hb outputs from rig A.
- **Decision P4-C: recompute them.** Storing a rig hash in `dataHistory` was ruled out by §1.1 ("dataHistory never references the ledger").
- **Implementation:**
  - `pipelineInfo` has a new field `rigDependent` (default `false` in `createPipelineInfo`; `true` for HemoCompute).
  - In `matchHistoryAndSetPreload`, a history match whose matched DAG nodes include a rig-dependent step is **not reused** when its `timestamp` is older than the folder's latest `acknowledgeConflict` event (`UMITRigLedger.lastAcknowledgementTime`, the ledger's UTC time converted to unzoned local).
  - Only two timestamps are compared; nothing is written to `dataHistory` or the ledger, so §1.1 holds.
  - Shorter matches that stop before the rig-dependent step are still reused. Downstream steps re-run because their input re-runs.
  - The ledger is read only when the pipeline contains a rig-dependent step. An unreadable ledger makes every rig-dependent output stale, so the step runs and its own guard reports the problem. An undated history entry also counts as stale.
  - Each refused reuse prints `[CACHE] Not reusing "<file>": it was computed before the folder's acknowledged rig conflict (...)`. The stale files stay on disk; new results get the usual `_N` suffix.
- **Test:** `testAcknowledgedConflictForcesRecomputationWhenSkipping` (HemoCompute → spatialGaussFilt chain):
  - an unchanged re-run reuses the branch;
  - after rebind + acknowledge, HbO is recomputed with different values;
  - the following run reuses the new outputs.
- **Pre-existing finding (unchanged):** a pipeline whose **leaf** is HemoCompute is never skipped at all, because HemoCompute's two outputs (HbO, HbR) make the full-branch match ambiguous (`numel(fileList) ~= 1`). Skipping only reuses HemoCompute when it sits inside a longer branch.

### Tests and validation
- `validateTaskChanges(ChangedFiles=[10 files], TestTargets=[test/Analysis/HemoCompute, TestCoregistrationRigResolution, TestUMITRigOptics, TestUMITRigLedger])` → **106/106 passed**.
- `test/PipelineManager` (subsystem, shared infrastructure) → **93/93 passed**.
- Static: `PipelineManager.m` has 74 advisories against 75 at HEAD (none new; the two `tfPreflight` ones introduced here were fixed). The other changed files have 0.

Brief Phase 4 test list:

| Required test | Covered by |
|---|---|
| Archived library rig still analyses; self-contained | `testArchivedLibraryRigStillAnalyses` (also deletes the rig folder and recomputes identically) |
| Library rig edited after commit has no effect | `testLibraryEditAfterCommitHasNoEffect` (identical HbO/HbR, no new events) |
| Modified filter set errors / new filter set allowed | `testModifiedFilterSetErrorsAndNewFilterSetIsAllowed` (also preflight, and no output rewritten); `testCommittedFilterSetMissingFromLibraryIsUsed` |
| Conflicted folder blocks; acknowledged works, hadConflict stays | `testConflictedFolderBlocksAndAcknowledgementUnblocks` |
| Errors for no ledger / NoRig | `testUnboundFoldersAreErrors`, `testUnboundDatasetIsRefusedNotDefaulted` |
| Commits only when not already committed | `testCommitsOpticsAndFilterSetOnce` |
| Batch: stop the affected folder, logged | `testPipelineSkipsOnlyTheFolderFailingRigPreflight` |

### Open items
- Manual: apply and check `PASTE_INTO_MLAPP.txt` (Coreg2Cams).
- `promptSelectLibraryRig` (the dialog) has no automated test; it is covered by the manual checks. Its choice logic (`buildLibraryRigChoices`) is tested.

## Phase 5 — DataViewer (backend done; app edits pending)

### Files
- `GUI/DataViewer/DatImageSource.m`: **write-on-read removed**. It no longer calls `ensureDatasetRigAssociation`; it only validates the folder ledger (read-only). `RigAssociationIssue` now reports an *invalid* ledger (warning `Umitoolbox:DatImageSource:invalidRigLedger`). A missing or "no rig" ledger is not an issue at load time.
- New, all read-only except where noted:
  - `GUI/DataViewer/describeFolderRigLedger.m`: display summary. States `Unknown`/`NoRig`/`Bound*`/`Conflicted`/`Invalid`; rig name; `hadConflict`; committed parts; label text and tooltip lines; action availability. Never throws for an invalid ledger.
  - `GUI/DataViewer/planFolderRigBinding.m`: Bind/Rebind preview, including the conflict preview **before** confirmation and the confirmation text. Archived rigs are refused (`Umitoolbox:DataViewer:archivedRig`).
  - `GUI/DataViewer/applyFolderRigBinding.m` (writes): re-plans, then `bind`/`rebind` with Origin `DataViewer`.
  - `GUI/DataViewer/formatRigLedgerHistory.m`: readable event history. Rigs are named as recorded in the ledger, so it reads the same on an empty library.
  - Dialogs: `runFolderRigBindingDialog.m` (picker with the Active rig pre-selected → preview → confirm → apply), `promptAcknowledgeRigConflict.m` (explains the conflict and what acknowledging means; optional comment), `showRigLedgerHistoryDialog.m`.
  - `UMITRigLedger.findSnapshot` (public, read-only).
- Paste 2 (`2_DataViewer.txt`): 1 Inspector change and 9 code edits; no new components, properties or callback bindings.
  - The "Assign Rig" submenu becomes "Folder Rig", still filled at run time, with Bind / Rebind Rig..., Acknowledge Rig Conflict... and Rig Ledger History... (enabled from the summary).
  - The folder label and tooltip come from the ledger summary, with no session-mirror fallback.
  - The load alert is now "Invalid Rig Ledger".
  - The dead methods `resolveRigDisplayName` and `getCurrentRigStatusText` are removed.
  - The coregistration button no longer passes the project session rig (C-6).
- Tests: new `test/DataViewer/TestFolderRigLedgerActions.m` (6 tests); `TestDatImageSource` rewritten for read-only loading (2 tests). `TestDataViewerManagementMenu` will expect 'Folder Rig' after paste 2 (a follow-up in the checklist).

### Design choice: DataViewer Bind/Rebind keeps the folder's coregistration binding
A Bind/Rebind from DataViewer takes the **optics** of the chosen library rig and keeps the folder's current **coregistration** snapshot (none for an Unknown/NoRig folder). The transform applied to the data can only change through a re-import (Phase 4, ImagesClassification). Taking the library rig's active transform here would bind a transform that was never applied, and would create a spurious coregistration conflict on dual-camera folders. Tested in `testRebindKeepsTheFolderCoregistration`. *To confirm in review.*

## Phase 6 — Cleanup and boundaries (project side done; app edit and deletions pending)

### Project side (minimal change; no redesign)
- `UMITProjectStore.addSession`: no longer calls `ensureDatasetRigAssociation` (which stamped the Active rig into `AcqInfos.mat`). The session's `rigUUID`/`rigID` are a **read-only mirror of the folder ledger's bound rig** (private `iFolderLedgerRig`), empty when the folder is not bound or the ledger is invalid. Nothing is written to the SaveFolder.
- Project validation: **`session_rig_required` removed**, since an empty session rig is valid (unbound folder). `session_rig_incomplete` (half-filled pointer) is kept.
- `migrateSessionRigAssignment`: now refreshes the mirror from the folder ledger (read-only towards the folder; no-op when unbound or already equal). It no longer stamps the Active rig.
- Error texts of `rigAssignmentNotSupported` (`addSession`, `updateSessionMetadata`) now point to DataViewer's Folder Rig menu.
- `GUI/promptAddSessionDialog.m` (used by ProjectBindingManager): the rig preview comes from the folder ledger (no Active-rig preview).
- Paste 3 (`3_ProjectManagerTool.txt`): the session rig display and the tool's own Add Session preview come from the ledger.
- Tests: `TestUMITProjectStore` fixture binds the folder ledger instead of `assignDatasetRig`; 3 tests rewritten (`testNewSessionMirrorsFolderRigLedger`, both `migrateSessionRigAssignment` tests). `TestProjectManagerTool` setup and `setupSyntheticRegistrationFolder` bind the ledger.

### Dead code removed now
- `UMITRigStore.ensureDatasetRigAssociation` (the write-on-read / Active-rig stamping) and its test. It had no caller left in `.m` files or in any `.mlapp` (all six apps grepped).
- Already removed in Phase 4: the HemoCompute Active-rig fallback.
- `assignDatasetRig` and `readDatasetRigAssociation` are marked **RETIRED** in their help. They are still called by DataViewer/ProjectManagerTool until pastes 2–3 are applied, and are deleted afterwards (see the checklist follow-ups). The same applies to `getOrCreateDefaultRig` (Coreg2Cams, paste 1).

### Definition-of-done audit (2026-10-14)
- Active-rig resolution in `.m` code: only `ImagesClassification` (dual-camera import; allowed) and `UMITRigStore` internals.
- In apps: DataViewer (5 hits), Coreg2Cams (1) and ProjectManagerTool (3) are exactly what pastes 1–3 remove. RigManagerTool, ProjectBindingManager and PipelineManagerTool have none.
- Read paths that wrote: `DatImageSource` (fixed). `addSession` and `migrateSessionRigAssignment` were project write paths that also wrote `AcqInfos.mat` (fixed).

### Validation (Phases 5–6 work)
- DataViewer/ledger: 68/68 (`TestFolderRigLedgerActions`, `TestDatImageSource`, `TestDatImageSourceReading`, `TestUMITRigLedger`).
- Subsystems: ProjectManagement + ImageRegistration + DataViewer + RigManagement → 321 run, 0 failed (4 incomplete: the known empty-library guards).
- After deleting `ensureDatasetRigAssociation`: `TestUMITRigStore` + ImageRegistration + ProjectManagement → 200 run, 0 failed (same 4 guards).
- Static: no new advisories (`UMITProjectStore` 30 vs 31 at HEAD; `DatImageSource` 6 vs 7; the new files 0).

## Manual validation round (2026-10-15)

### Fixes from Bruno's checks
- **Coreg2Cams check c (archived bound rig) failed:** `startupFcn` raised `uialert` while the figure was still hidden ("Figure handle 'Visible' value must be 'on'"). Paste `1b` moves the alert after `Visible = 'on'`; `promptSelectLibraryRig` no longer alerts on a hidden figure.
- **Component renamed:** Bruno renamed `AssignRigMenu` → `FolderRigMenu` (callback `FolderRigMenuSelected`) while applying paste 2. Kept; paste `2b` updated the method bodies. `TestDataViewerManagementMenu` and `MLAppDelta_tests` updated.
- **Rebind confirmation hidden:** the modal rig picker stayed open after OK. Its nested-function callbacks kept the workspace, and so the `onCleanup`, alive. `promptSelectLibraryRig` and `promptAcknowledgeRigConflict` now delete their dialog right after `uiwait`; a regression test was added and verified to fail without the fix.
- **Rig library hygiene:** archiving OiS200 for check c made a leftover test rig the Active rig (archiving the Active rig requires a replacement). With Bruno's approval OiS200 was reactivated and 14 leftover test rigs from interrupted runs (2026-10-02 to 2026-10-09) were deleted, after a verified zip backup.

### All app edits applied and promoted (2026-10-15)
- Pastes 1, 1b, 1c (Coreg2Cams), 2, 2b (DataViewer) and 3 (ProjectManagerTool) were applied and manually validated. All three references were promoted with `confirmAndUpdateReference`, and fresh deltas are empty.
- The ProjectManagerTool reference dated from 2026-08-03, so earlier changes that had been validated but never promoted were accepted with it (for example the rig dropdown replaced by a read-only display).
- **Deleted** now that no app or code calls them:
  - in `UMITRigStore`: `assignDatasetRig`, `readDatasetRigAssociation`, `getOrCreateDefaultRig`, and the private `iLoadAcqInfoStream`, `iSaveLoadedMatAtomic` and `iFindLegacyDefaultUUID`;
  - `GUI/resolveRigDisplayName.m` and `GUI/buildSessionRigChoices.m`;
  - their tests: 7 in `TestUMITRigStore`, 2 in `TestProjectManagerTool`.
- All 11 `.mlapp` files were scanned before deleting.
- Validation: RigManagement, ProjectManagement, DataViewer, DataImport, ImagesClassification, ImageRegistration and `MLAppDelta_tests` gave 364 run, 0 failed and 3 incomplete (the empty-library guards in `TestImportEmptyRigLibrary`).
- Static analysis: no messages.
- User docs and the error-message reference moved to their own task, `_tasks/active/task_docs_rig_ledger.md`. It is blocked until this redesign is accepted and Bruno's ProjectManager changes are done.

### Decision: no hashes in user-facing text (Bruno, 2026-10-15)
- Users are told *which rigs* are involved and that their optics or transform **differ**. Snapshot hashes are never shown, and there is no detailed list of differences.
- Covered:
  - the rebind conflict preview (`iDetectConflict`) and the conflicted-folder error;
  - `reimportConflict`, the coregistration and commit mismatch errors and `filterSetModified`;
  - the HemoCompute filter-set warning and the import console line;
  - DataViewer's bind confirmation, the folder-label tooltip and the Ledger History.
- History now says per rebind whether the optics and the coregistration were changed, unchanged, added or removed.
- Hashes remain only in the ledger file itself and in the ledger-integrity errors (`hashMismatch`), which are technical diagnostics.

### Decision: no opt-out of camera coregistration at import (Bruno, 2026-10-15)
- Bruno raised "dual-camera data cannot be imported without applying the calibration". An `ApplyCoregistration` import parameter was implemented, then **reverted on his request**: a wrong calibration is fixed with Recalibrate in DataViewer_Coreg2Cams, then a re-import.
- `ImagesClassification` is back to its Phase 3 behaviour. The guard in `TestRunImagesClassification/testPipelineInfo` still asserts there is no `ApplyCoregistration` parameter.

## Project-side touchpoints (left for a later project redesign)
Every place where project code touches rig data, with its state after Phase 6. These are not redesigned here: the session still keeps a rig pointer, now only as a mirror of the folder ledger.

| # | Touchpoint | State after Phase 6 | For the later redesign |
|---|---|---|---|
| C-1 | `getUMITProjectSchema.m`: session fields `rigUUID`/`rigID` (required fields, may be empty) | Kept; mirror of the folder ledger's bound rig | Drop the mirror, or make it derived on read instead of stored |
| C-2 | `UMITProjectStore` validation: `session_rig_required` / `session_rig_incomplete` | `session_rig_required` removed; `session_rig_incomplete` kept | Remove with C-1 |
| C-3 | `UMITProjectStore.addSession` → was `ensureDatasetRigAssociation` | Read-only mirror via private `iFolderLedgerRig`; no folder write | — |
| C-4 | `migrateSessionRigAssignment`, `iStampSessionRig` (its health-gate exemption still lists `session_rig_required`, harmless), `updateSessionMetadata` rig-field rejection | Mirror refresh from the ledger; rejection texts point to DataViewer | Mirror goes stale after a DataViewer rebind until `migrateSessionRigAssignment` runs (C-5) |
| C-5 | Mirror freshness: a DataViewer Bind/Rebind does not update project sessions | Unchanged (by design: the ledger is authoritative; the mirror is informational) | Derive on read, or refresh on bind |
| C-6 | DataViewer passed the session mirror RigUUID/RigID to Coreg2Cams | Removed in paste 2; Coreg2Cams resolves the folder ledger first anyway (paste 1) | — |
| C-7 | `promptAddSessionDialog.m`, `ProjectManagerTool.mlapp` (`createRigDisplay`, `promptSessionInfo`) rig display/preview | `.m` done; the app is paste 3 | — |
| — | `ProjectManagerTool.mlapp` `promptSessionInfo` duplicates `promptAddSessionDialog.m` | Unchanged (both now read the ledger) | Consolidate |
| — | Context structs (`context.rigUUID/rigID` in `resolveProjectBinding`/`findSessionByDataFolder`, UMITProjectStore ~L736, L840, L1013) | Unchanged; they carry the mirror | Remove with C-1 |
