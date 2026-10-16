# Implementation notes — project management redesign

Running log for the brief "Project Management redesign — optional projects, folder-committed image references, JSON project ledger" (2026-10-15). Context: `AUDIT_project_management.md`, `AUDIT_rig_management.md`, `IMPLEMENTATION_NOTES_rig_ledger.md`.

## Phase status

| Phase | Status |
|---|---|
| 0. Reconciliation and shared ledger core | Done, accepted (2026-10-15) |
| 1. Alignment correctness | Done; ImageAlignmentTool paste applied and checked (2026-10-16) |
| 2. Reference model | Done, awaiting review (backend; GUI items pending) |
| 3. Project ledger | Done (backend), validated; awaiting review |
| 4. Commitment, guards and fresh import | Done (backend), validated; accepted (2026-10-15) |
| 5. Remove project gates; GUI | Backend done, validated; GUI and docs pending (app edits) |
| 6. Store robustness | Done (backend), validated; GUI entry points pending |
| 7. Adoption | Done (backend and DataViewer offer, 2026-10-16), validated |

## Decisions confirmed with Bruno

| # | Question | Decision |
|---|---|---|
| §4 Q3 | `hadConflict` after a fresh import | **Keep it** (true by history). |
| P0-A | Existing rig ledgers holding the retired `reimport` event | **Invalid** (unknown event type → `invalidLedger`). No migration, pre-release. |
| P0-B | `assertCanReimport` / `reimportConflict` once a fresh import resets the optics commitment | **Removed.** A re-import with a different Active rig rebinds without conflict. |
| §4 Q1 | What may still run after apply | **Refuse both** estimate and apply, whatever the reference. Only a fresh import reopens alignment. |
| §4 Q4 | `alignment_bkp_<ts>` backup zip | **Removed.** |
| §4 Q5 | Session registration-transform resources | **Stop writing them.** `DataParams` is the only authority. The store API, `activeRegistrationTransformUUID` and the PMT handlers are removed in Phase 6 (with audit item 15). |
| P1-A | `UseFile='auto'` when the reference was not built from a `.dat` in the folder | Use the reference's `.dat` name if the folder has it; otherwise the folder's **only** Y-X-T `.dat`; several → error (`AutoSourceFileAmbiguous`), none → `AutoSourceFileNotFound`. Never guesses between channels. |
| §4 Q6 | `file` source kind | **UMIT exports + plain images** (TIF/MAT/UMT through `readImageReferenceSource`, kind `file`, no identity). |
| §4 Q7 | Adoption (Phase 7) | **Same hash in the library: show the match only** (nothing recorded). **No match: offer to add** the folder's reference to the subject library, with confirmation. `DataParams.registration` is never modified. |
| §4 Q8 | Existing `.umitlink` files | **Discard.** No converter; a folder with only a `.umitlink` reads as unbound (the file is ignored, never deleted automatically). |
| §4 Q9 | Dead PMT resource nodes | **Removed** (PMT paste, 2026-10-16). Also remove the registration-transform resource API (Q5). Reference store APIs stay. |
| P1-B | GUI work while Bruno cannot use App Designer (2026-10-15) | **Defer all `.mlapp` edits**; continue backend phases and collect GUI items under "Pending GUI work". |

| §4 Q2 | Fresh-import deletion scope | **Known UMIT artefacts only:** an explicit list of the patterns UMIT creates (survey at the end of this file). Unknown files stay and are listed in the confirmation. Raw input and both ledgers are always kept. |
| Q2-b | Backup before a fresh import deletes | **Keep the GENBACKUP / ERASE choice:** GENBACKUP zips everything that will be deleted; the zip is kept. |
| Q2-c | `importFromTif` into a non-empty folder | **Always a fresh import** (all channels of a recording are imported in one call). |
| P4-D | Stale project/ledger locks (2026-10-16) | **No menu item:** a "locked" error offers to clear the lock (`offerStaleLockClearing`; project lock >= 60 min, folder ledger lock >= 10 min; recent locks kept). Used in DataViewer, the Project Binding Manager, the Project Manager (`showException`) and the Folder Rig dialogs. A ledger "locked" error names its folder, which is used when it exists. |
| P4-C | One fresh-import decision (2026-10-16) | **`backupOpts` is the only zip-or-erase decision.** `AllowFreshImport` removed; `backupOpts` allows `GENBACKUP` (default) and `ERASE` only (no `REFUSE`). DataViewer's dialog only informs (lists the deletions; Continue/Cancel, `confirmFreshImport(plan, 'AskBackup', false)`). A direct call with `backupOpts = ''` keeps its zip/delete/cancel dialog on a desktop and means `GENBACKUP` in batch. **Deviation from brief 1.6:** pipelines are no longer refused for a used SaveFolder; the default `GENBACKUP` keeps the deleted files in a zip. `Umitoolbox:DataImport:saveFolderNotEmpty` and `PipelineManager.setParameterValue` were removed. |
| P6-A | Single RawFolder authority (audit item 7) | **Option A, derive** (2026-10-16): `DataParams.folders.RawFolder` is the authority; binding or relocating a RawFolder that differs is refused, a mismatch is a validation warning; `SessionInfo.rawDataFolder` stays as the record of where the RawFolder ledger is. |

All §4 questions are answered.

## Pending GUI work (deferred, 2026-10-15)
- **DataViewer — DONE (2026-10-16):** pastes 01–09 applied by Bruno, all manual checks passed, reference promoted (delta: 7 new, 14 modified, 3 removed methods), paste folder removed. Contents: project binding state from the ledger (`describeProjectBinding`, new state `stale`); orphan and stale-link recovery (Management > Project > Repair Folder Binding, Project Ledger History); binding problems reported in the Status Label on load; context label with `Align: aligned / not aligned`; gates P1/P2 removed for the Image Alignment Tool; fresh import on File > Import = an informative dialog, then the `backupOpts` parameter (P4-C); v2 reference preview; dead `startSetImageReference` removed; stale locks cleared only from a "locked" error (`GUI/DataViewer/offerStaleLockClearing.m`, also used by the Folder Rig dialogs). DataViewer items below are therefore done, except where noted.
- DataViewer `readExistingImageReference` must read `registration.reference.image` (v2). Until then, **Preview Dataset Registration Reference shows nothing for folders whose DataParams were rewritten as v2** (it reads `registration.referenceImage` defensively and finds none; nothing errors).
- DataViewer `startSetImageReference` (dead code per the audit) writes `registration.referenceImage/referenceFile/referenceDescription` through `updateDataParam`; with v2 those writes would fail validation. Delete it with the next DataViewer paste.
- **(DONE in DataViewer and ProjectBindingManager, 2026-10-16)** **DataViewer and ProjectBindingManager decided "bound" with `isfile(<binding file>)`** (`refreshProjectContextForCurrentData`, `inspectCurrentBinding`, `BindFolderButtonPushed`). With the ledger, a folder whose binding was removed (unbind, orphan removal) keeps its `UMITProjectLedger.json` (history), so until the app edit those apps show such a folder as **invalid** instead of unbound. Use `UMITProjectStore.getProjectBindingState(folder)` (`'unbound' | 'bound' | 'invalid'`). Bound folders and never-bound folders behave as before. Texts mentioning `.umitlink` in DataViewer (`classifyProjectBindingFailure` alerts) should be updated in the same paste.
- **DataViewer import menu (fresh import confirmation):** `PASTE_DataViewer` 04 + 06 (applied) + 07 (to apply): an informative dialog (`planFreshImport` → `confirmFreshImport(plan, 'AskBackup', false)`, Continue/Cancel), then the import parameter `backupOpts` decides zip or erase.
- **ImageAlignmentTool — DONE (2026-10-16), gate P3:** Select Image Reference offers Subject Library (bound folders only), UMIT .dat File and Reference File; every source goes through `resolveImageReferenceSource`, and Apply passes `'referenceSource'` (the recorded descriptor). The aligned-source warning is shown in an alert. Bruno's manual checks passed; reference promoted (delta: 3 modified methods).
- **Phase 5 gates:** P1, P2 (DataViewer) and P3 (ImageAlignmentTool) DONE 2026-10-16. **P5 DONE (2026-10-16):** Image Reference Manager standalone mode (create from the current frame or a file without a project, Export... as a reference file) and "Export as File..." for library references; DataViewer opens it in standalone mode when the folder has no project session. Manual checks passed; references promoted. **All Phase 5 gates (P1, P2, P3, P5) are removed.**
- **DataViewer status:** registration label/tooltip from `describeFolderRegistration`; project label from `describeFolderProjectLedger`; a "Project Ledger History" dialog from `formatProjectLedgerHistory` (like the rig one).
- **Binding states:** replace `classifyProjectBindingFailure` (message-text based) with `describeProjectBinding`; offer orphan recovery at load, stale-link clearing (`removeStaleProjectBinding`), `invalid`/`conflict` repair, and `clearStaleLock` (audit item 14).
- **ProjectBindingManager — DONE (2026-10-16):** pastes PBM and PBM-2 applied by Bruno, reference promoted (delta: 2 new, 4 modified methods). Bind/Rebind calls `bindFolderToProjectSession` (all or nothing, item 11); orphan rebind with a new session via `ReplaceOrphanBinding` (item 12); binding state from the ledger (`getProjectBindingState`, new private `currentBindingState`); a "locked" error offers `offerStaleLockClearing` (result in the status line). **Decision PBM-1 (Bruno):** the dialog binds to a NEW session only, in Assign and Rebind mode; existing sessions are no longer listed (binding to one moved that session off its own SaveFolder). Moving an existing session to another folder is done in the Project Manager; the backend still allows it (the store's `sessionAlreadyBound` guard, item 13). Manual checks a–d by Bruno; e–g run by Claude through the app's callbacks in an isolated `UMIT_ROOT` with copies of MultiChannel_SmallData: orphan → session whose folder is gone (bound, rebind event), orphan → healthily bound session (refused, nothing changed), bound-then-unbound folder opens as unbound, stale ledger lock cleared then bind succeeds, recent lock kept.
- **ProjectManagerTool — DONE (2026-10-16):** paste PMT applied by Bruno, reference promoted (delta: 9 modified methods; 10 methods, 4 callbacks and 4 buttons removed). Dead resource code removed (item 15, §4 Q9): `addResourceNodes`, the resource editor case, `saveResourceChanges`, Set Active / Archive / Restore / Open Location handlers and buttons (`onSetActiveResource` was the last caller of `setActiveRegistrationTransform`), `buildResourceFileFilter`, `promptResourceInfo`, `resourceNodeText`, and `onRetrySaveFolder` (no button). Subject resource registries still show as read-only tables; Add Image Reference stays. Sessions show a "Folder Rig" row from the SaveFolder's rig ledger (T8; a legacy `rigID` is ignored). `showException` offers `offerStaleLockClearing` on any "locked" error, which now reads the locked folder from the ledger error itself (new test `testFolderIsTakenFromTheLockError`). The Add Session rig tooltip no longer says sessions store a rig (also in `GUI/promptAddSessionDialog.m`). Manual checks a–c by Bruno; d–f run by Claude through the app's callbacks in an isolated `UMIT_ROOT` with copies of MultiChannel_SmallData: Add Session tooltip; stale project lock on Save cleared then saved, recent lock kept; stale folder-ledger lock on Remove Session named the folder, was cleared, and the session stayed until the retry removed it.
- **Adoption offer — DONE (2026-10-16):** paste ADOPT applied by Bruno, DataViewer reference promoted (delta: 1 new, 4 modified methods). New `GUI/DataViewer/offerReferenceAdoption.m` (test `test/DataViewer/TestOfferReferenceAdoption.m`): plan → match shown only (`inLibrary`), or a dialog Add to Library / Add and Make Active / Not Now (`canAdd`); a lock error offers lock clearing. DataViewer runs it quietly at the end of `launchProjectBindingManager` only when the binding changed (new private `folderBindingKey`: projectUUID/sessionUUID from the ledger), and from the new menu item Management > Project > Add Folder Reference to Library... (enabled when the binding is valid). Manual checks a–d by Bruno (on `MultiChannel_SmallData_AdoptionCheck`, a copy aligned to the original's green.dat without a project); e–g run by Claude in an isolated `UMIT_ROOT` (unaligned folder: no offer after binding, menu says "not aligned"; "Use Binding" on a bound aligned folder: no offer, confirmed by a positive control; menu disabled with no data or an unbound folder).
- **Documentation** (audit §11.6 pages, error-message reference): projects optional, reference sources, commitment, fresh import, project ledger; written with the app edits.

## Phase 0 — Reconciliation and shared ledger core

### 0.1 Rig/project reconciliation report (current state, 2026-10-15)

What rig Phase 6 already changed on the project side (verified in code):
- `UMITProjectStore.addSession` no longer writes `AcqInfos.mat`: the session's `rigUUID`/`rigID` are copied read-only from the folder rig ledger by the private `iFolderLedgerRig` (`UMITProjectStore.m:1237`, `:6573-6592`), empty when the folder is unbound or the ledger invalid. (Audit R3/W2 and item 22 are therefore already fixed.)
- `session_rig_required` was removed from validation; `session_rig_incomplete` (half-filled pointer) remains (`:5685-5701`).
- `migrateSessionRigAssignment` (`:1976-2003`) refreshes the mirror from the ledger through `iStampSessionRig` (`:4799-4840`); it no longer stamps the Active rig. Its health-gate exemption still lists the now-nonexistent code `session_rig_required` (`:4818`, harmless).
- `ensureDatasetRigAssociation`, `assignDatasetRig`, `readDatasetRigAssociation`, `getOrCreateDefaultRig`, `GUI/resolveRigDisplayName.m` and `GUI/buildSessionRigChoices.m` are deleted (audit R10 and the duplicated display helper are gone).
- `GUI/promptAddSessionDialog.m` and `ProjectManagerTool.mlapp` (`createRigDisplay`) display the rig from `describeFolderRigLedger` (ledger), not from the session mirror.
- DataViewer no longer passes the session rig to Coreg2Cams; its rig label and Folder Rig menu come from `describeFolderRigLedger` (audit R7, R8 resolved).

Remaining touchpoints and inconsistencies (to resolve in Phase 6, §1.12):

| # | Touchpoint | Location | Resolution planned |
|---|---|---|---|
| T1 | Session schema fields `rigUUID`/`rigID` (required, may be empty) | `getUMITProjectSchema.m:252-253`, comment `:289-292`, header `:25`, `:184` | Remove (R1) |
| T2 | `session_rig_incomplete` validation | `UMITProjectStore.m:5685-5701` | Remove with T1 |
| T3 | `addSession` rig stamp (mirror) and rejection of rig fields | `:1138-1141`, `:1172`, `:1237`, `:1281-1282` | Remove the stamp; keep rejecting unknown fields generically (R3) |
| T4 | `updateSessionMetadata` rig-field rejection text pointing at `migrateSessionRigAssignment` | `:1898-1923` | Simplify once the fields are gone (R4) |
| T5 | `migrateSessionRigAssignment` + `iStampSessionRig` | `:1976-2003`, `:4799-4840` | Delete (R5) |
| T6 | Context structs carry the mirror (`context.rigUUID/rigID`) | `resolveProjectBinding` `:736-737`; `:840-841`; `findSessionByDataFolder` `:1013-1014` | Delete; no remaining consumer in `.m` code or apps (DataViewer uses `describeFolderRigLedger`) |
| T7 | Store help text mentions the rig pointer | `UMITProjectStore.m:19` | Update |
| T8 | PMT displays a "Rig" row only because session metadata has a `rigID` field (`ProjectManagerTool.mlapp` › metadata grid `if strcmp(fieldName,'rigID')`, `createRigDisplay`; session `preferred` field order; resource table preferred names) | app | Strict-mode paste in Phase 6: render the ledger rig row explicitly instead of keying it on the schema field |
| T9 | `TestUMITProjectStore` mirror tests (`testNewSessionMirrorsFolderRigLedger`, two `migrateSessionRigAssignment` tests) | tests | Rewrite/delete with T1–T5 |
| T10 | `DataParams.cameraCoregistration.rigID/resourceUUID` | `createDataParams.m:199-207`, `validateDataParams.m:447-504`, `loadDataParams.m:134` | Not project-owned (rig provenance, best-effort, written by ImagesClassification). Out of scope here; listed for completeness |
| T11 | `UMITProjectStore` → `UMITRigLedger.validate` dependency (only through `iFolderLedgerRig`) | `:6584` | Disappears with T3 |

No other inconsistency was found between the rig implementation and project code: `PipelineManager`, `EventsManager`, `Analysis/**` (except `createRegistrationTform`) and the import functions have no project dependency, and the rig code has none on project code.

### 0.2 Shared ledger core

**File:** `subFunc/UMITLedgerCore.m` (static-method class, no state). Placed in `subFunc` (shared I/O helpers) so the import functions and both ledgers depend on it without depending on `UMITProjectStore` or `UMITRigStore`.

**Spec.** Every I/O, lock and input-check function takes a ledger spec from `UMITLedgerCore.makeSpec(struct(...))`: `fileName`, `lockFolderName`, `format`, `schemaVersion`, `errorPrefix` (e.g. `Umitoolbox:RigLedger`), `label` (e.g. `rig ledger`), `clearLockCall`, optional `bodyFields` (top-level fields between `createdOn` and `events`), `folderArgName`, `folderLabel`. This keeps each ledger's error identifiers and message texts.

**API**
- Hashing: `canonicalJSON(value, hashOptions)`, `computeHash(value, hashOptions)` (`hashOptions.listFields`, `hashOptions.errorPrefix`), `hashText`, `hashBytes` (SHA-256 of raw bytes, `sha256:` prefix; for the Phase 2 content hash), `shortHash`, `quantize`.
- Events and time: `newLedger(spec, body)`, `newEvent(type, origin)`, `appendEvent`, `newUUID`, `nowText`, `dateText`, `localTime`.
- File: `filePath`, `decodeFile` (`<prefix>:readFailed`), `normalizeHeader` (required top-level fields, events as a cell list), `cellOfStructs`, `headerIssues` (format, schemaVersion, ledgerUUID), `eventIssues` (non-empty, base fields, unique `eventUUID`), `writeLedger(folder, spec, ledger, verifyFcn)` (temp file → caller verification → backup → rename → restore on failure; `<prefix>:writeFailed`).
- Locking: `acquireLock(folder, spec, operation)` (atomic `java.io.File.mkdir`; `lockInfo.mat` with operation, createdOn, processID; `<prefix>:ledgerLocked` / `lockCreationFailed`), `clearStaleLock(folder, spec, minimumAge)` (`<prefix>:lockNotStale`), `lockAge`.
- Messages and inputs: `buildErrorMessage`, `folderText`, `existingFolder`, `requireOrigin` (`<prefix>:invalidInput`), `isText`, `maybeText`, `nullable`, `deleteIfPresent`, `removeFolder`.
- Constants: `HashPrefix`, `SignificantDigits` (6), `ShortHashDigits` (12), `ToolboxVersion` ('2.0.0', now the single copy; `UMITRigLedger.ToolboxVersion` points to it), `TimestampFormat`, `BaseEventFields`.
- Core-own error IDs: `Umitoolbox:LedgerCore:invalidSpec`, `Umitoolbox:LedgerCore:invalidInput` (`hashBytes`), and `Umitoolbox:LedgerCore:hashFailed` when no `errorPrefix` is given.

**Canonical form.** Unchanged from the rig ledger (normative text in the class help and in `IMPLEMENTATION_NOTES_rig_ledger.md`), generalised: the fields always written as lists are a caller option (`listFields`; the rig passes `cameras`, `illuminations`). The exact canonical-text test passes unchanged on both classes.

**`UMITRigLedger` refactor (no behaviour change).** Its public API, file format, error identifiers and message texts are unchanged; the generic private helpers now delegate to the core. Two internal differences, both invisible to callers:
- The duplicate-`eventUUID` check moved from replay into the integrity check (`eventIssues`); both run before any state is derived, the message text is the same, and the outcome is still `invalidLedger`.
- Lists of missing fields in integrity messages are in schema order instead of alphabetical order.

### 0.3 `freshImport` replaces `reimport` (rig ledger)
- `UMITRigLedger.recordFreshImport(saveFolder, 'Origin', origin)` writes a `freshImport` event: `previousOpticsHash`, `previousCoregistrationHash`, `previousFilterSets` (list of `{id, hash}`; informational, not validated). Replay resets the optics, coregistration and filter-set commitments. The binding (or NoRig) is kept, so the folder becomes BoundUncommitted (or stays NoRig) and an open conflict ends. Earlier events stay; `hadConflict` stays true (Q3). A folder without a ledger → `noLedger`.
- `recordReimport`, `assertCanReimport` and `reimportConflict` are removed (P0-B). A ledger with a `reimport` event is invalid (P0-A).
- `ImagesClassification`: on a re-import it writes `freshImport` before classifying; the dual-camera preflight no longer refuses different optics. Help text updated. **Interim until Phase 4:** only the rig side of a fresh import exists. Old outputs are removed only as far as `genBackupFolder` already does on a re-import; artefact deletion, confirmation, the pipeline refusal, `importFromTif` and the project-ledger event come in Phase 4.
- `GUI/DataViewer/formatRigLedgerHistory.m`: history line for `freshImport`. No `.mlapp` references the removed names (all 11 apps' `document.xml` grepped).
- `_tasks/active/task_docs_rig_ledger.md`: the import-page bullet updated to the fresh-import behaviour.

### Tests and validation
- New `test/subFunc/TestUMITLedgerCore.m` (15 tests): exact canonical text (incl. null, matrix rows, control-character escape), list fields, key-order independence and 6-digit stability, JSON round trip, SHA-256 known vectors (`abc`, empty), caller error prefixes, new ledger/event shape (no user/host), header and event integrity failures (foreign format, version, missing UUID, empty events, missing base fields, duplicate UUIDs), malformed/missing files, atomic write with no leftovers, rejected candidate keeps the previous file, lock contention and release, stale-lock clearing, input checks.
- `TestUMITRigLedger`: the re-import section replaced by 6 fresh-import tests (reset of all commitments with the binding kept and history retained; conflict ended with `hadConflict` kept; NoRig stays NoRig; legacy `reimport` invalid; drop-coregistration after a fresh import; no ledger → `noLedger`). `testAssertCanReimportChecksOnlyCommittedOptics` and `testReimportResetsOnlyCoregistration` removed with the behaviour they covered.
- `TestRunImagesClassificationCoregistration`: event lists use `freshImport`; `testReimportWithDifferentOpticsAfterCommitIsRefused` replaced by `testFreshImportWithDifferentOpticsAfterCommitRebinds`. `TestRunImagesClassification`: single-camera re-import keeps the binding and is BoundUncommitted; `testRepeatedSingleCameraImportRecordsFreshImport`.
- Before the `freshImport` change, the refactor alone: `TestUMITRigLedger` 33/33.
- After all changes, `validateTaskChanges` (subsystem) over the core tests, `test/RigManagement`, `test/Analysis/{DataImport, ImagesClassification, importFromTif, HemoCompute}` and the DataViewer ledger tests: **257 run, 0 failed, 3 incomplete** (the known empty-library guards of `TestImportEmptyRigLibrary`, which only run on an empty rig library). Static: 3 pre-existing FNDSB/LOGSUM/ISCL advisories in `ImagesClassification.m`, untouched code; 0 elsewhere.
- `test/ProjectManagement`, `test/Analysis/ImageRegistration`, `test/PipelineManager`: **239/239 passed**.

### Open issues (Phase 0)
- The empty-library import guards were not re-run on an empty library (requires moving the real rig library aside; same procedure as rig Phase 2). The change does not touch the no-rig paths.
- `ToolboxVersion` is still a manual mirror of `Umit.prj` (now in one place).

## Phase 1 — Alignment correctness

### Files
- `subFunc/applyImageAlignmentToFolder.m`: the canonical apply backend, rewritten around a staged commit.
- `subFunc/assertFolderCanBeAligned.m` (new): the shared guard used by estimate and apply, GUI and pipeline.
- `Analysis/ImageRegistration/applyRegistrationTformOnFolder.m`: now a pipeline wrapper around the backend.
- `Analysis/ImageRegistration/createRegistrationTform.m`: guard, apply record untouched, resize composed into the transform, `UseFile='auto'` fix.
- `PASTE_INTO_MLAPP.txt`: 3 strict-mode edits to `ImageAlignmentTool.mlapp` (applied 2026-10-16; file removed).
- Tests: `TestApplyRegistrationTformOnFolder`, `TestCreateRegistrationTform`, `TestImageAlignmentCharacterization`, `TestSameSizeRewritersCharacterization`.

### One apply path (§1.7)
- `applyImageAlignmentToFolder` is the only apply implementation. `applyRegistrationTformOnFolder` loads `DataParams`, runs the guard, shows its confirmation dialog and calls the backend with the stored `tform` and `referenceImage` and `keepEstimate=true`.
- `keepEstimate`: the pipeline applies exactly the transform `createRegistrationTform` stored, so the estimate fields (method, reference, source, QC, managed provenance) are kept and only the apply record is written. The GUI (`keepEstimate=false`) rebuilds the estimate fields from its own inputs and clears stale QC fields. Nothing is carried over from an earlier estimate in either path, including `imageReferenceUUID`/`Checksum` and `resourceUUID`.
- Both paths transform `.dat`, image `.umt`, current `.roi`, `DataParams.mask` and `DataParams.view`, and write `manualAlignmentReport_<ts>.mat`. Tested: `testGuiAndPipelineApplyGiveIdenticalFolderState` (data, ROI vertices and masks, `DataParams` mask and view are identical).

### Staged commit replaces the backup zip (§4 Q4)
The zip was also the GUI path's failure rollback, so removing it needed a replacement:
1. Everything is validated first: DataParams, `.dat` layouts and sizes, `.umt` and `.roi` sizes, the view size, the precomputed mask, and the new `DataParams` (`validateDataParams`).
2. **Stage:** each transformed file and the report are written to `<name>.aligning.tmp`. The originals are untouched, so a failure here only deletes the temp files.
3. **Commit:** each original is moved to `<name>.aligning.bak` and its staged copy installed, then `saveDataParams` (validated, atomic) runs. Any failure restores every original from `.bak`, removes newly created files and deletes the temp files. After success the `.bak` files are deleted.
4. Either way, a failure writes `ImageAlignmentFailedReport_<ts>.txt` with its stage.
- The raw `save` of `DataParams` (audit §6.3 `:226`) is gone: the backend loads with `loadDataParams` and saves with `saveDataParams`.
- Disk use: all transformed files exist side by side during the commit (as the zip did).
- Tested: `testFailedCommitRestoresEveryOriginal`. On Windows, `red.dat` is held open so its move fails mid-commit; every file is byte-identical afterwards, no temp or backup file remains, and the folder can still be aligned.

### Guard and commitment (§4 Q1, audit item 17)
- `assertFolderCanBeAligned` → `Umitoolbox:ImageRegistration:alreadyAligned` when `registration.isRegistered`. It is called by the backend (so the GUI path is guarded too; before, it could re-align an aligned folder) and by `createRegistrationTform` before anything is resolved or written. The message names the operation, the folder, the reference (description), when and by what it was applied, and the fix (fresh import). Phase 4 adds the source kind and short hash.
- `createRegistrationTform` no longer writes `isRegistered`, `isReviewed`, `applied*`, `confirmationMode` or `resourceUUID`.
- `AllowReapply` was removed from `applyRegistrationTformOnFolder`; passing it is now an `inputParser` error.

### Behaviour changes in the pipeline apply (consequence of §1.7, not separately asked)
- **Interpolation: `nearest` → `linear`** (the GUI backend's, which the brief makes canonical). Pipeline results differ numerically from before. The frozen `applyRegistrationTformOnFolder` case in `TestSameSizeRewritersCharacterization` was retired; the backend stays pinned by `TestImageAlignmentCharacterization`, whose recorded values are unchanged.
- Output frame = the reference image size (before: the input size, which had to equal it).
- ROIs, mask and view are now transformed.
- The pipeline's `.dat` layout support (any Y-X-leading layout, NaN restoration) moved into the backend. For the GUI this adds Y-X, Y-X-E, Y-X-F and Y-X-T-E files, and NaN pixels no longer spread through linear interpolation. Data without NaN give the same values as before.
- Error IDs of refused inputs are now the backend's: `ImageAlignmentTool:UnsupportedLayout`, `:UmtSizeMismatch`, `:ViewSizeMismatch` (new), `:RoiSizeMismatch` (new), `:CommitFailed` (new). The wrapper keeps `MissingDataParams`, `MissingTform`, `UnknownReferenceSize` and `Cancelled`.

### Latent GUI-path bugs found and fixed
- **`.umt` files written by `saveData` were skipped.** `saveData` saves with `-struct` (the variables are the UMT fields), but the backend only looked for a named UMT variable and reported "no UMT structure found". Both layouts are now handled; the `-struct` layout is rewritten in the same format as `save2umt` (`-v7.3`).
- **Current `.roi` files were skipped as "legacy".** The backend checked an older ROI layout (`fileType`, top-level `imageSizeYX`). It now recognises `getROISchema` files (`schemaName 'UMIT_ROI'`). It updates `imageInfo.imageSizeYX`, `view.originXY_px` and `pixelSize_px_per_mm`, clears `statsImage.image` and the intensity stats, recomputes the geometric stats, and validates with `validateROIFile` before staging.
- **Image `.umt` entries of another frame size** were warped silently; now `UmtSizeMismatch`.
- New consistency check: `DataParams.view.imageSizeYX`, when set, must equal the `.dat` frame size (`ViewSizeMismatch`).

### Transform convention (`createRegistrationTform`)
- The stored transform now always maps the folder's own image coordinates into reference coordinates. If the moving image had to be resized to the reference for the estimate (only possible with an empty `view.imageSizeYX`), the resize `u -> s(u - 0.5) + 0.5` is composed back in (stored as `affinetform2d`). A different aspect ratio cannot be a similarity transform → `Umitoolbox:createRegistrationTform:AnisotropicResize`. QC metrics still describe the estimated transform.
- Both paths therefore share the backend's size convention, and the wrapper's old reference-size gate is no longer needed. Tested: `testDifferentSizeMovingImageIsComposedIntoTransform` (128×128 recording onto a 64×64 reference) and `testDifferentAspectRatioIsRefused`.

### `UseFile='auto'` (audit item 20, P1-A)
- Removed IDs: `AutoSourceFileUnavailable`, `AutoSourceFileNotDat`. New: `AutoSourceFileAmbiguous`. `AutoSourceFileNotFound` now means "no Y-X-T `.dat` in the folder".
- Tested with TIF, MAT and UMT reference sources (`testAutoWithTifMatOrUmtReferenceUsesTheOnlyYXTDat`) and the never-guess rule (`testAutoWithoutUsableSourceNeverGuessesBetweenChannels`).

### ImageAlignmentTool (§4 Q5, audit item 16) — applied 2026-10-16
`PASTE_INTO_MLAPP.txt` removes `registerRegistrationTransformResource` and its call, so no transform resource is written, before or after apply. It also stops passing `resourceUUID` and updates the confirmation and progress texts (no backup; aligned only once). Until the paste is applied, the backend still accepts `resourceUUID` and writes it as given (never carried over). The option and the `DataParams` field go away with `DataParams.registration` v2 (Phase 2).

### Tests and validation
- `validateTaskChanges` (targeted) on `test/Analysis/ImageRegistration` and `test/subFunc/rewriters`: **86/86**, 0 static messages. Pre-existing Code Analyzer messages in `applyImageAlignmentToFolder.m`, and one in `TestSameSizeRewritersCharacterization.m`, were cleaned so the changed files are message-free.
- Subsystem run (ImageRegistration, rewriters, DataParams, PipelineManager, DataViewer, ProjectManagement, Analysis/PipelineManager): **353/353**.

| Required test (brief Phase 1) | Covered by |
|---|---|
| Identical folder state from GUI and pipeline apply | `testGuiAndPipelineApplyGiveIdenticalFolderState` |
| Atomic write on failure | `testFailedCommitRestoresEveryOriginal`, plus the existing preflight no-mutation tests |
| Re-estimate refused on an applied folder | `testEstimateIsRefusedOnAnAlignedFolder`; re-apply: `testAlignedFolderRefusesEveryReapply`, `alignedFolderRefusesASecondAlignment` |
| `UseFile='auto'` with TIF, MAT and UMT references | `testAutoWithTifMatOrUmtReferenceUsesTheOnlyYXTDat` |

### Open items (Phase 1)
- Manual (2026-10-16): Bruno applied the paste; the delta showed exactly the 3 edits; checks a–d passed; the reference was promoted (`confirmAndUpdateReference`).
- User docs for `createRegistrationTform`, `applyRegistrationTformOnFolder` and the manual registration tool (no `AllowReapply`, no backup zip, aligned once, linear interpolation) → Phase 5 documentation pass.
- `ImageAlignmentTool.getAlignmentImpactSummary` / `confirmRoiImpactBeforeApply` were not reviewed against the ROI-schema fix; paste check d covers the visible behaviour.

## Phase 2 — Reference model (backend)

### Files
- New in `subFunc/`: `computeImageContentHash`, `resolveImageReferenceSource`, `referenceSourceDescriptor`, `validateReferenceSourceDescriptor`, `newRegistrationRecord`, `upgradeRegistrationRecord`, `validateRegistrationRecord`, `describeRegistrationReference`, `verifyRegistrationReference`, `exportImageReferenceFile`, `importImageReferenceFile`, `clearImageReferenceIdentity`.
- New in `GUI/DataViewer/`: `normalizeImageReferenceStruct`.
- Changed: `genImageReferenceStruct`, `validateImageReferenceStruct`, `readImageReferenceSource` (non-interactive options), `UMITProjectStore` (import check), `createDataParams`, `loadDataParams`, `saveDataParams`, `validateDataParams`, `updateDataParam`, `applyImageAlignmentToFolder`, `applyRegistrationTformOnFolder`, `createRegistrationTform`, `assertFolderCanBeAligned`.
- Tests: new `test/ImageReference/TestImageReferenceResolver.m` (8 tests), new store test `testImageReferenceImportChecksIdentityAndVersion`; updated `TestDataParams`, `TestCreateRegistrationTform`, `TestApplyRegistrationTformOnFolder` and the `TestUMITProjectStore` reference fixture (now identity-less).

### Content hash (normative, §1.3)
`computeImageContentHash(image)`:
1. Normalize as `genImageReferenceStruct` does: `single`, non-finite → 0. Also −0 → +0, so the bytes do not depend on the sign of zero.
2. Hash input, in order: the UTF-8 bytes of `single` followed by one `0x00` byte; the number of dimensions as `uint32` little-endian; each dimension of `size(image)` as `uint64` little-endian; the image values as `single` little-endian, column-major. Big-endian hosts swap bytes first.
3. `'sha256:' + lowercase hex` via `UMITLedgerCore.hashBytes`.

The same image gives the same hash from any source and after any re-save (tested: a project library copy vs a `.dat` extraction; two saves of one reference have different file checksums but the same hash). The hash is computed when a reference is selected or created. Estimate and apply re-verify the stored image (`verifyRegistrationReference` → `Umitoolbox:ImageRegistration:referenceHashMismatch`; the apply backend calls it in pipeline mode, and Phase 4 extends it to the estimate path). `validateDataParams` only checks the hash format, so loading stays cheap.

### `ImageReference` v2 (§1.8)
- `version = 2`, new field `contentHash` (validated against the image).
- Identity is optional. `projectUUID`/`projectName` and `subjectUUID`/`subjectID` are both-or-neither pairs, and a subject needs a project. `genImageReferenceStruct` now needs only `Name`.
- Version 1 is still accepted by `validateImageReferenceStruct` (with its original strict identity) and upgraded in memory by `normalizeImageReferenceStruct`, so existing libraries keep working.
- **Store import (audit item 21):** `addImageReference` accepts only v2 (`Umitoolbox:UMITProjectStore:imageReferenceVersion`). Identity fields present in the payload must name this project and the owner subject (`...:imageReferenceIdentityMismatch`). The check runs in `iAddManagedResource` after the owner is resolved, so a failed check writes nothing.
- Image Reference Manager already builds payloads through `genImageReferenceStruct`, so it produces v2 without an app change.

### Reference resolver (§1.4)
`resolveImageReferenceSource(source, 'Interactive', tf)` returns `image`, `contentHash`, `name`, `description`, `axisConvention`, `imageSizeYX`, `source` (descriptor) and `sourceFile` (the file the image came from, for `UseFile='auto'`). It reads only: it never writes to the source folder, the target or the store (tested with folder and store checksums).

| Kind | Input | Descriptor fields filled |
|---|---|---|
| `project` | `store`, `subjectID`, optional `resourceUUID` (default: the active reference) | `path` (store file at selection time), `projectUUID`, `projectName`, `subjectUUID`, `subjectID`, `resourceUUID`. The registry checksum is verified, and the payload identity must match the library that holds it. |
| `umitDat` | `datPath`, or `saveFolder` + `datFile` | `path` (.dat), `saveFolder`, `extractionRule='evenlySpacedMean'`, `extractionRuleVersion=1`, `sourceAligned`, `sourceCommittedHash`, `sourceCommittedName`, `warning` |
| `file` | `path`, optional `variable`/`entry` | `path`, `fileFormat` (`umitReference` for an exported file, `image` for TIF/MAT/UMT), identity fields of an exported file, `entry` (MAT variable / UMT entry used) |

- **Descriptor:** one flat struct for every kind (`referenceSourceDescriptor`); fields that do not apply stay empty. Kind `legacy` exists only for records upgraded from v1.
- **Extraction rule (`umitDat`):** the mean of up to 100 evenly spaced frames, exactly as `readImageReferenceSource` reads a `.dat` (unchanged code).
- **Validity (`umitDat`):** the `.dat` folder must hold an `AcqInfos.mat` that loads and contains `AcqInfoStream`; otherwise `Umitoolbox:ImageRegistration:invalidDatSource`. A `.dat` given as kind `file` is refused (use `umitDat`).
- **Aligned source:** allowed with the warning `Umitoolbox:ImageRegistration:alignedSource`. The text says the source data are already transformed, that the target will be interpolated a second time, and gives the source's committed reference (name + short hash). It is also stored in `source.warning`, so the folder record keeps it. The source folder is read with `loadDataParams` (read-only).
- **Non-interactive mode:** `readImageReferenceSource` gained `Variable`, `Entry` and `Interactive`. With `Interactive=false`, several candidates raise `readImageReferenceSource:AmbiguousSource` instead of opening `listdlg`. Pipelines (and `createRegistrationTform`) use it.
- Error IDs: `Umitoolbox:ImageRegistration:invalidSource`, `:invalidDatSource`, `:referenceNotFound` (also wraps store lookup errors), `:referenceChecksumMismatch`, `:referenceIdentityMismatch`, `:alignedSource` (warning).

### `DataParams.registration` v2 (section 2)
```
registration.version    = 2
registration.reference  = [] | struct(image, contentHash, name, description, axisConvention,
                                      imageSizeYX, selectedOn, source)
registration.estimate   = [] | struct(tform, transformType, method, qcStatus, qcWarning, qcMetrics,
                                      qcFigureFile, qcPreviewImageFile, createdOn, movingSource,
                                      referenceHash)
registration.apply      = struct(isRegistered, isReviewed, appliedOn, appliedBy, confirmationMode,
                                 appliedReferenceHash, manualAdjustment, notes)
```
- **Removed:** `referenceFile` (→ `reference.source.path`), `imageReferenceUUID` and `imageReferenceChecksum` (→ `reference.source.resourceUUID` + `reference.contentHash`), `resourceUUID` (§4 Q5), and the flat fields (`referenceImage` → `reference.image`, `tform` → `estimate.tform`, and so on). `manualAdjustment` is now in the schema (`apply.manualAdjustment`).
- **Rules (`validateRegistrationRecord`):**
  - each sub-record has exactly its fields;
  - hash formats are checked;
  - `reference.imageSizeYX` must equal the image size;
  - `isReviewed` requires an estimate;
  - when not aligned, the apply fields stay empty;
  - when aligned, an estimate exists, `appliedReferenceHash` equals `reference.contentHash` and the reference size equals `view.imageSizeYX`.
  - `estimate.referenceHash` may differ from the current reference: re-selection is free before apply, and the Phase 4 apply guard handles it.
- **Version 1 on disk (deviation: in-memory upgrade).** The brief requires no migration, but every existing save folder has a v1 `DataParams.mat` and v2-only validation would make all of them unloadable. `loadDataParams`, `saveDataParams` and `createDataParams` therefore upgrade a v1 record in memory (`upgradeRegistrationRecord`):
  - `referenceImage` → `reference` (hash computed, name from the old description, source kind `legacy` with the old `referenceFile` as path and `imageReferenceUUID` as `resourceUUID`);
  - `tform` + QC → `estimate`;
  - apply fields → `apply` (`appliedReferenceHash` = the reference hash when aligned).

  Loading never writes (tested). The file becomes v2 the next time a writer saves it, and a legacy estimate can still be applied (tested).
- **Writers:**
  - `createRegistrationTform` writes `reference` (from the resolver, kind `project`) + `estimate` and never touches `apply`.
  - The apply backend writes `apply` and commits `appliedReferenceHash = estimate.referenceHash`. In pipeline mode (`keepEstimate`) it refuses an estimate made against another reference (`ImageAlignmentTool:EstimateReferenceMismatch`). In GUI mode it rebuilds `reference` and `estimate` from the GUI inputs.
  - New backend option `referenceSource` (descriptor). By default it is built from `imageReference` (kind `project` when the struct carries a project identity, else `file`).
  - The options `imageReferenceUUID`/`imageReferenceChecksum` from Phase 1 are gone; `resourceUUID` is accepted and ignored until the ImageAlignmentTool paste is applied.
- `updateDataParam` no longer has the v1 `registration.tform` → `isRegistered` side effect.
- `assertFolderCanBeAligned` now names the committed reference with its source kind and short hash (`describeRegistrationReference`). Hashes appear in alignment messages because the brief asks for them (§1.5); this differs from the rig rule of no hashes in user-facing text.

### Export and import of standalone references (§1.8)
- `exportImageReferenceFile(ImageReference, path, 'ClearIdentity', tf)`: writes a MAT file with the single variable `ImageReference` (v2), atomically through `saveMatAtomic`. It works with or without a project, and the file is a `file` source.
- `importImageReferenceFile(store, subjectID, path, 'ClearIdentity', tf, 'ResourceInfo', s)`: imports into a subject library through `addImageReference`, so the version and identity checks apply. `ClearIdentity` adopts a reference made in another project or subject. The original file is never modified.
- Creating a reference without a project from any source (GUI project-less mode) is `genImageReferenceStruct(reference.image, 'Name', ...)` + export; the Image Reference Manager mode itself is Phase 5 (GUI).

### Tests and validation
Run in chunks (one combined run crashed MATLAB; every chunk passed separately, so the crash was the run's size, not a test):
- DataParams, ImageReference (including the resolver), ImageRegistration, rewriters, ProjectManagement: **201/201**;
- subFunc, RigManagement and the resolver: **130/130**;
- DataImport, ImagesClassification, importFromTif, HemoCompute: **116 passed, 3 incomplete** (the known empty-library guards);
- PipelineManager (both folders): **109/109**;
- DataViewer: **60/60**.

Static: 0 messages in the new and changed files; the count in `UMITProjectStore.m` is unchanged from HEAD (30, all pre-existing).

| Required test (brief Phase 2) | Covered by |
|---|---|
| Each source kind | `testUmitDatSourceIsReadOnlyAndDescribed`, `testFileSourcePlainImages`, `testExportedReferenceWithAndWithoutProject`, `testProjectAndDatSourcesGiveTheSameHash` |
| An invalid `.dat` folder | `testUmitDatNeedsAReadableAcqInfos` |
| The aligned-source warning | `testAlignedSourceIsAllowedWithAWarning` |
| Source folder byte-identical after selection | `testUmitDatSourceIsReadOnlyAndDescribed`, `testAlignedSourceIsAllowedWithAWarning` |
| Same image from a project and a `.dat` → same hash | `testProjectAndDatSourcesGiveTheSameHash` |
| Hash stable across re-saves | `testHashIsContentBasedAndStableAcrossResaves` |

### Open items (Phase 2)
- GUI items listed under "Pending GUI work".
- `createRegistrationTform` still takes its reference only from the bound project; the source parameter and the resolution order come in Phase 4, with the selection writer and the apply guard.

## Phase 3 — Project ledger (backend)

### Files
- New `ProjectManagement/UMITProjectLedger.m` (static class on `UMITLedgerCore`; it depends on neither `UMITProjectStore` nor `UMITRigStore`, so import code can use it).
- `ProjectManagement/UMITProjectStore.m`: every binding read and write now goes through the ledger (below).
- `ProjectManagement/getUMITProjectSchema.m`: `files.projectBinding = 'UMITProjectLedger.json'`; `metadataVariables.projectBinding` removed. The in-memory binding record (`projectBinding.version`, `requiredProjectBindingFields`) is unchanged.
- `Analysis/funcTemplateAcquisitionInitializer.m`: `UMITProjectLedger.json` counts as dataset content (as `.umitlink` did).
- `Analysis/ImageRegistration/createRegistrationTform.m`: help text.
- Tests: new `test/ProjectManagement/TestUMITProjectLedger.m` (9 tests). `TestUMITProjectStore`: 23 binding-existence assertions moved from `isfile(.umitlink)` to `getProjectBindingState`; the reciprocity and orphan tests check the ledger; new `testProjectLedgersOfRawAndSaveFolders`.

### Format (`UMITProjectLedger.json`)
- Header from the shared core: `format = 'umit-project-ledger'`, `schemaVersion = 1`, `ledgerUUID`, `createdOn`, `events`. No body fields.
- Atomic writes (temp file, verification by replay, rename) and the folder lock `UMITProjectLedger.lock` (`ledgerLocked` / `lockNotStale`), as in the rig ledger.
- **Binding** = exactly `{bindingUUID, projectUUID, subjectUUID, sessionUUID, folderRole}`: lowercase canonical UUIDs, `folderRole` ∈ {`processedDataFolder`, `rawDataFolder`}. These are the only fields UMIT reads and validates (strict: unknown or missing fields → invalid).
- **Events** (base fields from the core):

| Type | Fields | Allowed when |
|---|---|---|
| `bind` | `binding`, `snapshot` | not bound |
| `rebind` | `binding`, `previousBinding`, `machine`, `snapshot` | bound; `previousBinding` must equal the current binding |
| `unbind` | `previousBinding`, `reason` | bound |
| `orphanRemoval` | `previousBinding`, `machine` | bound |
| `freshImport` | — | bound |

- **Snapshot** (`makeSnapshot(ProjectInfo, subjectID, sessionID)`): `takenOn`, `project` (every `ProjectInfo` field except `subjectRegistry`, datetimes as ISO-8601 UTC text), `subjectID`, `sessionID`. It is frozen and never read back, validated or refreshed; any content is tolerated, including retyped or added fields (tested).
- **States** (derived): `Unknown` (no file), `Bound`, `Unbound` (the binding was removed; the file and its history stay). Orphan, mismatch and conflict are still derived by `resolveProjectBinding` from the binding plus the store.
- `machine` is `COMPUTERNAME` (or `HOSTNAME`, or the Java host name). It is recorded only where the brief asks for it (rebind; also orphan removal, which records which machine found the project absent). As in the rig ledger, no user name is recorded.
- Error identifiers: `Umitoolbox:ProjectLedger:invalidLedger`, `alreadyBound`, `notBound`, `invalidBinding`, `ledgerLocked`, `lockNotStale`, `lockCreationFailed`, `readFailed`, `writeFailed`, `invalidInput`.

### Deviation: an `unbind` event (not in the brief's event list)
The brief lists bind, rebind, orphan removal and freshImport. Several store operations used to **delete** the `.umitlink` of a folder: removing a session or subject, deleting a project, relocating a session away from a folder, and unbinding a RawFolder. Deleting the ledger would erase its history ("history is never rewritten"), so these record an `unbind` event with a `reason`: `sessionRemoved`, `removedFromProject`, `relocated` or `folderUnbound`.

### Store changes (persistence only; the rollback logic is unchanged)
- `readProjectBinding` reads the ledger and returns the same binding record as before (`createdOn` / `modifiedOn` come from the event times). A folder that is not bound raises `bindingReadFailed` as a missing link did.
- New public `getProjectBindingState(folder)` → `'unbound' | 'bound' | 'invalid'`, for the apps (pending GUI work).
- Installing a binding (`iInstallBinding`) writes `bind` for an unbound folder and `rebind` for a bound one: an orphan being replaced, `moveSessionToSubject`, `rebindSessionToProject`, a repair or relocation onto a bound folder, or a co-located folder that stays bound for the other role. A rebind to an identical binding writes nothing, so a co-located raw and save folder that share one ledger get one event, not two.
- Removing a binding (`iRemoveBinding`) writes `unbind`. An invalid ledger cannot take an event: it is left as it is, with the warning `Umitoolbox:UMITProjectStore:invalidLedgerKept` (the old code deleted a malformed link).
- `removeOrphanProjectBinding` writes `orphanRemoval` (it no longer deletes anything).
- **Rollback:** each operation backs the ledger bytes up first (`iBackupLedger`) and restores them if it fails (`iRestoreLedger`, or the existing recovery copies via `copyfile`). The failed operation never happened, so this is a rollback, not a history rewrite. Restores became unconditional, because an unbound ledger file still exists.
- "Has a binding" checks (`bindingFileExists`, the validation codes `missing_project_binding` / `invalid_project_binding`, `findSessionByDataFolder`, runtime status) use the ledger state: `iHasBinding` is true for a bound folder or an invalid ledger, which both block a new binding as an existing link did.
- §4 Q8: an old `UMITProjectBinding.umitlink` is ignored (the folder reads as unbound) and never deleted (tested).
- Code Analyzer: `UMITProjectStore.m` keeps its 30 pre-existing messages (same count as HEAD).

### Tests and validation
- `test/ProjectManagement` (with the new ledger tests): **100/100**.
- Before the new tests: ProjectManagement + ImageRegistration + ImageReference **167/167**.
- The chunk that hung the MCP desktop MATLAB (DataViewer + DataParams + ImageRegistration + the ledger test) was re-run in fresh **batch** MATLAB processes (`matlab -batch`, per-run timeout). All of it passes there (batches C and E below). The hang was not reproduced, so the blocking test was not identified. A desktop-only modal dialog is the likely cause; `TestDataViewerManagementMenu` passes alone in batch (2/2).

| Required test (brief Phase 3) | Covered by |
|---|---|
| Bind, rebind and orphan-removal histories | `testBindRebindUnbindHistory`, `testOrphanRemovalIsAnEvent`, `TestUMITProjectStore/testDetectAndRemoveOrphanProjectBinding` |
| Snapshot fields ignored by validation | `testSnapshotContentIsNeverValidated` |
| Identity mismatch detected | `testIdentityIsStrictlyValidated` (malformed identity), `testProjectLedgersOfRawAndSaveFolders` (valid identity that does not match the store → `resolveProjectBinding` error and `invalid_project_binding`) |
| RawFolder ledger when Raw ≠ Save | `testProjectLedgersOfRawAndSaveFolders` |
| Reading never writes | `testReadingNeverWrites`, `testProjectLedgersOfRawAndSaveFolders` |

## Q2 survey — what a SaveFolder holds and what re-imports do today (2026-10-15, read-only)

**Files UMIT creates in a SaveFolder:**

| Category | Patterns |
|---|---|
| Raw input (never deleted) | `*.bin`, `*info.txt`, `Comments.txt`, `Snapshot*.png`, TIFF input `*.tif/*.tiff`, `info.json` |
| Import output | channel `*.dat`, `AcqInfos.mat`, `ROI.mat` (sub-ROI crop), `tformDualCam.mat`, legacy sidecars `<base>.mat`/`<base>_info.mat` |
| Folder state | `DataParams.mat` |
| Derived | analysis `*.dat`/`*.umt`; `*.roi`; `events.mat`; `dataHistory.mat`; `pipeLog.mat`; `pipeErrorLog_*.csv`; DataViewer ROI exports `*_ROI_*.csv`; `registrationQC_*.fig/.png`; `manualAlignmentReport_*.mat`; `ImageAlignmentFailedReport_*.txt`; `*MotionCorrection.mat` / `*MotionCorrectionQC.png`; `std_speckle.tiff`; ConvertToTiff exports `<base>.tif`, `<base>_C*_R*.tif`, `<base>_info.txt`; genBackupFolder zips `bkp_*.zip` / `<name>.zip` |
| Temp leftovers | `PMTMP_*`, `*_writing.*`, `*_compute.dat`, `*.dat.tmp`, `*.aligning.tmp/.bak`, `tp*.mat`, `tp*.json` |
| Ledgers | `UMITRigLedger.json` (+ `.lock/`), `UMITProjectLedger.json` (+ `.lock/`), legacy `UMITProjectBinding.umitlink` |
| Legacy names | `Text_events.mat`, `LogBook.mat`, `ImagingReferenceFrame.mat` |

**Re-import behaviour today:**
- **`ImagesClassification`:** writes the rig `freshImport` event, then `genBackupFolder(..., 'eraseFolder', true)`. That manages only top-level `*.dat`, `*.mat`, `*.roi`, `*.umt`: they are zipped (GENBACKUP) or deleted (ERASE). Everything else stays: QC figures, `.txt`/`.csv`, zips, temp leftovers, ledgers. Its "fixed" raw list never overlaps the managed patterns.
  - Derived exports named like raw files (`*.tif`, `*_info.txt`, `std_speckle.tiff`) would survive any pattern-based rule.
  - **Bug:** if the interactive backup dialog is cancelled, nothing is deleted but the import still overwrites the channels, and the rig `freshImport` event has already been written.
- **`importFromTif`:** no backup, no deletion, no ledger event. It overwrites or adds the channels it imports (old channels stay listed), keeps `DataParams` (so stale alignment survives), and asserts that the raw size, rate and binning match an existing `AcqInfos.mat`.
- **PipelineManager:** a folder counts as fresh only when it has no `AcqInfos.mat`. A non-empty folder with a valid `AcqInfos.mat` is not refused. The run likely re-saves `pipeLog`/`dataHistory` loaded before the import step, so history the import just deleted may be recreated (inferred from the call order, not run).
- **DataViewer:** no prompt before importing into a non-empty folder; `.bin` imports use GENBACKUP/ERASE silently.

## Phase 4 — Commitment, guards and fresh import (backend)

### Files
- New in `subFunc/`: `planFreshImport`, `applyFreshImport`, `prepareFreshImport`, `confirmFreshImport`, `selectRegistrationReference`, `normalizeReferenceSourceInput`.
- Changed:
  - `ImagesClassification`, `importFromTif`: fresh import, new `AllowFreshImport` parameter, and `backupOpts` (new in `importFromTif`);
  - `createRegistrationTform`: `ReferenceSource`;
  - `applyImageAlignmentToFolder`: apply guard;
  - `applyRegistrationTformOnFolder`, `assertFolderCanBeAligned`: requested reference in the message;
  - `PipelineManager`: reloads the history after an importer step.
- Tests:
  - new `test/Analysis/DataImport/TestFreshImport.m` (6 tests);
  - `TestRunImagesClassification`: new `testFreshImportResetsTheFolder`, and re-imports now pass `AllowFreshImport`;
  - `TestRunImagesClassificationCoregistration`: re-imports pass `AllowFreshImport`;
  - `TestCreateRegistrationTform`: 5 new tests.

### Fresh import (§1.6, §4 Q2)
- **Plan** (`planFreshImport`, read-only): a folder holds a previous import when it has a deletable artefact or a rig or project ledger.
  - **Deleted** (top level): only the known artefacts listed in the Q2 survey, plus every output recorded in `dataHistory.mat` (field `file`) and temporary leftovers.
  - **Kept:** raw files, both ledgers and their locks, a legacy `.umitlink`, earlier `*.zip` backups, subfolders and unknown files. The unknown ones are listed so a confirmation can show them.
  - ConvertToTiff exports and legacy sidecars are recognised by the folder's own channel names, so raw TIFF inputs are not confused with exports.
- **Apply** (`applyFreshImport`):
  1. GENBACKUP / name: zip the files to be deleted (`bkp_<ts>.zip`, made unique; kept). ERASE: no zip.
  2. Ledger events: rig `freshImport` when the rig ledger exists (all commitments reset); project `freshImport` only when the folder is bound. Unbound folders get no project ledger.
  3. Delete the files.
- **Consent** (`prepareFreshImport`), decided before anything is touched:
  - (Superseded 2026-10-16, decision P4-C: `AllowFreshImport` and the pipeline refusal were removed; `backupOpts` GENBACKUP (default) / ERASE decides.) Originally: `AllowFreshImport = true`: proceed (`backupOpts ''` means GENBACKUP).
  - Direct call with `backupOpts ''` on a desktop: `confirmFreshImport` lists the deletions and the kept unknown files and offers *Zip, delete and import* / *Delete and import* / *Cancel*. Cancel → `Umitoolbox:DataImport:freshImportCancelled`.
  - (Superseded 2026-10-16, P4-C: no refusal any more.) Originally: otherwise (pipelines, batch MATLAB, explicit `backupOpts` without the option): `Umitoolbox:DataImport:saveFolderNotEmpty`.
- **Importers:**
  - Both call `prepareFreshImport` after their read-only preflight (invalid rig ledger; dual-camera rig resolution) and before writing.
  - `ImagesClassification` no longer calls `genBackupFolder`. This **fixes the survey bug**: a cancelled dialog used to skip the deletion but import anyway.
  - The calibration path (`ApplyRig = false`) resets the folder without ledger events.
  - `importFromTif` gained `backupOpts` and `AllowFreshImport`. Any import into a folder holding an import is a fresh import (Q2-c), so the old "existing AcqInfos must match" rule no longer applies to re-imports.
  - Both pipelineInfo entries declare `backupOpts` (`GENBACKUP` default, `ERASE`); `AllowFreshImport` was removed on 2026-10-16 (P4-C).
- **Deviation (non-empty definition):** commit `7077569` let importers ignore stale leftovers when `AcqInfos.mat` is missing. A folder that holds only stale known artefacts (for example `.dat` files without `AcqInfos.mat`) now counts as a previous import: pipelines re-import it with `backupOpts` (zip first by default), and the confirmation deletes those leftovers. PipelineManager's own "fresh SaveFolder" rule (no `AcqInfos.mat`) is unchanged.
- **PipelineManager:** after a step whose function is an `acquisition-initializer`, the run reloads `dataHistory.mat` from disk (pruning entries of deleted files). The history loaded before the import is then not written back over the fresh import. No prompts: the importers refuse rather than ask, and guard errors stop only their folder (the existing per-folder catch).

### Reference source and guards (§1.4, §1.5)
- `createRegistrationTform` `'ReferenceSource'`: a source struct, or a path (`.dat` → `umitDat`, else `file`; `normalizeReferenceSourceInput`), which pipelines can pass as a parameter. Resolution order:
  1. the parameter;
  2. else, for a bound folder (`getProjectBindingState`), the subject's active reference;
  3. else `Umitoolbox:ImageRegistration:noReferenceSource`, which explains both ways to supply one.

  Projects are optional: an unbound folder estimates and applies with a `.dat` source (tested).
- `selectRegistrationReference(SaveFolder, source)`: records the selected reference (validated, atomic save), for the GUI.
  - Before apply it is free; an existing estimate is kept but then mismatches.
  - After apply, the same hash is a no-op and any other is `Umitoolbox:ImageRegistration:referenceCommitted`.
- **Apply guard** (shared backend, before any write): the estimate must have been made against the selected reference's hash, else `Umitoolbox:ImageRegistration:estimateReferenceMismatch` (both hashes named; fix: estimate again, or select the estimated reference). A missing estimate gives `...:noEstimate`.
- **After apply** (§4 Q1): estimate and apply are refused whatever the reference (`alreadyAligned`). The message names the operation, the folder, the committed reference (name, source kind, short hash), the requested one (same or different), and the fix: a fresh import. Guards compare hashes only: project, `.dat` and file sources follow identical rules (tested both ways).
- The aligned-source warning stays the only warning-level path. In pipelines it is a MATLAB warning (logged by PipelineManager's run output) and processing continues.

### Error identifiers added in Phase 4
`Umitoolbox:DataImport:freshImportCancelled` (`saveFolderNotEmpty` removed 2026-10-16), `:freshImportBackupFailed`, `:invalidInput`; `Umitoolbox:ImageRegistration:noReferenceSource`, `:referenceCommitted`, `:estimateReferenceMismatch`, `:noEstimate`.

### Brief Phase 4 test list
| Required test | Covered by |
|---|---|
| Estimate against A, then B, apply B (allowed) | `testEstimateAgainstAThenBThenApplyB` |
| Apply after the reference changed since estimation (error) | `testApplyAfterTheReferenceChangedIsRefused` (folder unchanged) |
| Post-apply attempt with a different reference (error; message content) | `testDifferentReferenceAfterApplyIsAnError` |
| Fresh import removes all listed artefacts, resets rig commitments, appends both events, allows a new alignment | `testFreshImportResetsTheFolder`, `TestFreshImport` |
| Fresh import of an unbound folder writes no project ledger | `testRepeatedSingleCameraImportRecordsFreshImport`, `testGenBackupZipsThenDeletesAndLogsTheRigLedger` |
| Pipeline import into a non-empty folder without the option (refused) | `testFreshImportResetsTheFolder` (refusal step), `testPipelineWithoutConsentIsRefusedBeforeAnyChange` |
| A batch from one `.dat` source shares one hash | `testBatchFromOneDatSourceSharesOneHash` |
| Project- and `.dat`-sourced folders follow identical rules | `testDifferentReferenceAfterApplyIsAnError` (`.dat` committed, project requested), existing project-sourced guard tests |

### Phase 3–4 validation (2026-10-15, batch MATLAB, after all Phase 4 changes)
The MCP desktop MATLAB was unresponsive, so validation ran in separate `matlab -batch` processes through `validateTaskChanges` (subsystem mode, per-run timeout), using a scratchpad runner (`runValidatorBatch.m` / `runWithTimeout.ps1`).

| Batch | Targets | Result |
|---|---|---|
| C | ImageRegistration, ImageReference, DataParams, ProjectManagement, TestFreshImport, subFunc (rewriters, ledger core, datHeader, …) | 213 run, 212 passed, 1 failed. The failure was `testMissingBindingDoesNotMutateFolder`, whose expectation was outdated: an unbound folder without `ReferenceSource` now raises `noReferenceSource`. Fixed and re-run in D2. |
| D | DataImport, ImagesClassification, importFromTif, RigManagement, HemoCompute | 233 run: 11 tests of the old re-import behaviour failed; 3 incomplete are the empty-library guards (on this machine the rig library has 4 rigs). The failing tests were updated, see below. |
| D2 | ImagesClassification, importFromTif, DataImport, TestCreateRegistrationTform (re-run) | 94 run, **91 passed, 0 failed**, 3 incomplete (the empty-library guards) |
| E | PipelineManager (both folders), DataViewer | **169/169** |

**Static analysis:** no new messages. `UMITProjectStore.m` (30), `PipelineManager.m` (74) and `ImagesClassification.m` (3) keep their pre-existing counts. Unneeded pragmas in the new or changed tests were removed.

**Tests changed for the new fresh-import rules:**
- The `testBackupErase…` / `testBackupGenBackup…` tests in `TestImagesClassification` and `TestRunImagesClassification` now pass `AllowFreshImport`. The ERASE ones also check the refusal without consent: stale leftovers make the folder a previous import.
- `TestImportFromTif`: 7 tests of the old additive or compare-AcqInfos behaviour were replaced by `testSecondImportIsAFreshImport`, `testMalformedPreviousAcqInfosIsReplaced` and `testImportIntoRigBoundFolderKeepsTheBinding`. The "existing AcqInfos.mat must match" code in `importFromTif` is now unreachable (a previous import is always deleted first) and can be removed in a cleanup.
- Refusal checks always pass an explicit `backupOpts`; an empty one would open the confirmation dialog in a desktop MATLAB and block the run.

**Store hygiene:** after the stopped run, a read-only listing of the real stores showed no leftovers from today. Older test leftovers remain (five "Image Registration Test …" projects from 2026-10-02 and rigs `CoregTargetPicker_4d0dbff0`, `CoregTargetPicker_6b2d9f5c`). They were not touched; the test store-root override is Phase 6 (audit item 24).

## Phase 5 — Remove project gates; GUI (backend only)

All GUI items of this phase are `.mlapp` edits, deferred (P1-B); they are listed under "Pending GUI work".

### Gates (audit §6.2)
| Gate | Where | Status |
|---|---|---|
| P1 `ensureProjectBinding` before alignment | DataViewer.mlapp | Pending app edit. The backend no longer needs a project. |
| P2 project-context builders | DataViewer.mlapp | Pending app edit. |
| P3 "No project context" in reference selection | ImageAlignmentTool.mlapp | Pending app edit (source selection: `resolveImageReferenceSource`, `selectRegistrationReference`). |
| P5 Image Reference Manager needs a store | ImageReferenceManager.mlapp | Pending app edit (project-less mode: `genImageReferenceStruct` without identity + `exportImageReferenceFile`, both ready). |
| P6 identity fields required in the reference struct | `genImageReferenceStruct` / `validateImageReferenceStruct` | **Removed** in Phase 2 (identity optional, v2). |
| P8 `createRegistrationTform` only from the project | `createRegistrationTform` | **Removed** in Phase 4 (`ReferenceSource`). |
| P7 managed library needs a store | `UMITProjectStore` | Stays, as the brief says. |

### New backend for the DataViewer status (`GUI/DataViewer/`, read-only, never throw for bad files)
- `describeFolderRegistration(saveFolder)`: alignment status from `DataParams.registration`. States `NoDataParams`, `Invalid`, `NoReference`, `ReferenceSelected`, `Estimated`, `EstimateStale` (estimated against another reference), `Aligned`. Gives the reference name, source kind, source path (human aid), short hash (the committed one when aligned), the aligned-source warning recorded at selection, apply time and user, and a display text plus tooltip lines. When the folder is not bound to a project, the tooltip adds the note that projects let one reference be reused across recordings.
- `describeFolderProjectLedger(folder)`: ledger state (`Unknown`, `Bound`, `Unbound`, `Invalid`), the bound identity and the names recorded in the latest snapshot. The names are display-only: the snapshot is never validated, and any content is tolerated (tested).
- `formatProjectLedgerHistory(folder)`: one line per event with the snapshot names of each bind/rebind, the machine of rebinds and orphan removals, and the reason of unbinds.

### Scenarios 2 and 2b without a project store
`test/Analysis/ImageRegistration/TestProjectlessAlignmentScenarios.m` runs on a UMIT root of its own (tag `noProjectStore`):
- scenario 2: one-off `.dat` source, estimate and apply, then a refused re-estimate;
- scenario 2b: a batch of three folders from one `.dat` source share one hash;
- an unbound folder without a source is refused.

Every test checks that the `projects` folder never appears below the root (the store is not even resolved) and that the source folder stays byte-identical. Scenario 1 (project source) stays covered by `TestCreateRegistrationTform` and `TestApplyRegistrationTformOnFolder`. `setupSyntheticRegistrationFolder` gained `'WithStore', false` (plain imported dataset, no rig or project).

### Not done in Phase 5
- User documentation (audit §11.6 pages) and the error-message reference: they describe the GUI flows, so they are written with the app edits. Listed under "Pending GUI work".

## Phase 6 — Store robustness

### Test store-root override (audit item 24; also resolves DFR-20261005-001's root problem)
- `getUmitFolder` honours the environment variable `UMIT_ROOT` (absolute folder; empty = not set; relative → `Umitoolbox:UmitFolder:invalidRootOverride`). It is the root of the project store, the rig library and the preferences. The help text now gives the real default path (`<USER>/Documents/LabeoTech/umIT`).
- `test/testHelpers/UMITTemporaryRootFixture.m`: a shared test fixture that points `UMIT_ROOT` at a fresh temporary folder and restores it afterwards. **Every** test class (107) declares it (`SharedTestFixtures`); `test/TestUmitRootIsolation.m` fails if a class does not. Classes with the same tag share one root while they run consecutively; `TestImportEmptyRigLibrary` (tag `emptyRigLibrary`) and `TestProjectlessAlignmentScenarios` (tag `noProjectStore`) get roots of their own.
- Effect: the empty-rig-library guards now run instead of being skipped (0 incomplete), and a listing of the real `Documents/LabeoTech/umIT` (paths, sizes, modification times) is identical before and after every validation batch of this phase.
- Left as is: the "snapshot and restore the Active rig" helpers still work (now on the temporary root). Deleting the old leftovers in the real stores needs your confirmation (see the report).

### Per-session validity and cheaper validation (§1.11, audit items 8, 11 of §6.2)
- Every validation issue has a new field `storeLevel`: false when it lies inside a registered subject's folder (subject or session errors), true otherwise (project metadata, schema, required folders, registries, folders below `subjects` that belong to no registered subject, duplicate UUIDs).
- `open` makes a project **read-only only for store-level errors**. A broken session no longer blocks the rest of the project.
- Mutations validate **what they touch**: `validate` gained a `'Scope'` option (subjects, sessions, subject trees, by UUID). Gates run in quick mode (no checksum of every resource) on the project level plus the touched nodes, before and after the write; the global UUID scan is skipped. Any error found blocks the mutation: `Umitoolbox:UMITProjectStore:invalidProject` when store-level (the store becomes read-only), otherwise the new `Umitoolbox:UMITProjectStore:invalidTarget` (the project stays writable). Scopes: `addSubject`, `updateProjectMetadata` → project only; subject edits → the subject node; session edits, bind/unbind, relocation → the session and its subject node; `removeSubject`, `renameSubjectID` → the subject with all its sessions; `moveSessionToSubject` → the session and both subjects; resources → their owner; `rebindSessionToProject` → the session in the source and the target subject.
- The narrow repair path (`relocateSaveFolder`, `repairSaveFolderBinding`, `removeSessionFromProject`) keeps its allowed codes on the selected session; `missing_acq_infos` became `not_a_dataset`.
- A scoped validation does not replace `LastValidationReport` (the GUI keeps the full report from `open`/`validate`).
- Resource lookup by UUID and the session-path search skip unreadable subject/session metadata instead of failing for the whole project.
- Tested: `testBrokenSessionBlocksOnlyItself`, `testMutationsValidateOnlyWhatTheyTouch`, `testUnregisteredFolderOpensReadOnly` (store-level stays read-only).

### Rig and project reconciliation (§1.12; Phase 0 list T1–T11)
| # | Resolution |
|---|---|
| T1 | `rigUUID`/`rigID` removed from `requiredSessionFields` and from new sessions. Existing session files that still hold them stay valid (extra fields are tolerated). |
| T2 | `session_rig_incomplete` removed. |
| T3 | `addSession` writes no rig. Unknown fields, rig fields included, are refused generically (`Umitoolbox:UMITProjectStore:unsupportedField`). |
| T4 | `updateSessionMetadata`: the rig-specific block is gone; non-editable fields are refused by the generic rule (`updateSessionFailed`). |
| T5 | `migrateSessionRigAssignment` and `iStampSessionRig` deleted. |
| T6 | Context structs (`resolveProjectBinding`, `getSessionContextByUUID`, `findSessionByDataFolder`) carry no rig fields. |
| T7 | Store and schema help texts updated. |
| T8 | PMT "Folder Rig" row from the SaveFolder's rig ledger: DONE (PMT paste, 2026-10-16). |
| T9 | Mirror tests removed or rewritten (`testSessionsKeepNoRig`). |
| T10 | `DataParams.cameraCoregistration.rigID/resourceUUID`: rig provenance, not project-owned; unchanged. |
| T11 | The store no longer depends on `UMITRigLedger`. |

### Registration-transform resources retired (§4 Q5/Q9)
- Removed: `addRegistrationTransform`, `setActiveRegistrationTransform`, `getActiveRegistrationTransform`, the `registrationTransform` resource type, `activeRegistrationTransformUUID`, and the session `transforms/registration/{active,archive}` folders for new sessions.
- Existing data stays valid: `schema.retiredResourceTypes = {'registrationTransform'}`. A record of a retired type is a **warning** (`retired_resource_type`), never listed or used; the session `transforms` folder is an allowed legacy child whose content is not checked (tested with a hand-made legacy session).
- `ImageAlignmentTool` no longer calls `addRegistrationTransform` (Phase 1 paste applied); no app calls the retired transform API (checked against the promoted references, 2026-10-16).

### Single RawFolder authority (audit item 7) — chosen: derive
- `DataParams.folders.RawFolder` (written by the importers) is the authority. `bindRawDataFolder(subject, session, '')` takes the RawFolder from it; a different folder is refused (`Umitoolbox:UMITProjectStore:rawFolderMismatch`), and so is a relocation. `SessionInfo.rawDataFolder` stays, as the record of where the store put the RawFolder's ledger.
- Validation warns (`raw_folder_mismatch`, not an error) when the two differ. A SaveFolder whose DataParams has no RawFolder (`'Missing'`, older imports) accepts any RawFolder.
- **Confirmed by Bruno (2026-10-16): option A.** Alternatives considered: B (SessionInfo central for bound sessions, syncing DataParams on bind/relocate) and C (two unreconciled records). SessionInfo.rawDataFolder is not removed: the store needs it to know where it wrote the RawFolder ledger, and unbound folders have no SessionInfo.

### Bindings: states, stale links, PBM back end (audit items 9, 11–14)
- `ProjectManagement/describeProjectBinding(folder)` (read-only): `unbound`, `invalid`, `orphan`, `conflict`, `projectUnreadable`, `stale`, `pathMismatch`, `bound`, with the fitting repair actions. It decides from the ledger, the store's metadata and error identifiers, **never from message text** (item 14). `isOrphanProjectBinding` no longer parses the message either (`iFindProjectRoots`), and the message-text helper `iIsProjectNotFoundError` is gone.
- **Stale links (item 9):** a link to an available project that no session matches any more (session removed, relocated or rebound while the folder was offline) can now be cleared with `UMITProjectStore.removeStaleProjectBinding(folder)`. It writes an `unbind` event with reason `staleBinding` (history kept); any other state → `Umitoolbox:UMITProjectStore:bindingNotStale`.
- **Item 11 (PBM transactional):** `ProjectManagement/bindFolderToProjectSession(folder, selection)` creates the project and/or subject when asked, then adds or binds the session. On failure it removes what it created (subject, then project), so a retry never finds a half-created project and never creates a duplicate. PBM calls it (done 2026-10-16).
- **Item 12:** `addSession(..., 'ReplaceOrphanBinding', true)` accepts a folder bound to a project that is not on this machine (a `rebind` event).
- **Item 13:** relocating a session onto an orphan folder (`ReplaceOrphanBinding`) is refused while the session's current folder still binds it (`Umitoolbox:UMITProjectStore:sessionAlreadyBound`); it is allowed once that folder is gone.
- Item 14's GUI entry points (orphan recovery at load, `invalid`/`conflict` repair, `clearStaleLock`) are app edits; the backend for all of them now exists.
- Item 15 (dead PMT resource nodes): removed (PMT paste, 2026-10-16).

### `deleteProject` and listings (audit item 23)
- `deleteProject` works on corrupted or offline projects, as documented: unreadable subject/session metadata, offline folders and folders whose ledger no longer binds them are **skipped and reported** (`result.skipped`, folder + reason) instead of blocking. A skipped folder keeps its ledger, which then reads as an orphan (clearable with `removeOrphanProjectBinding`).
- `resolveProjectRoot` and `listProjects` skip every internal folder whose name starts with `.` (`.__creating_*`, `.deleteProjectStage_*`, `.deleteProjectRecovery_*`).

### Smaller items (audit item 24)
- **IDs:** a new or renamed subject/session ID that is not a valid name is **refused** (`invalidID`) with the valid form in the message ("Use "Mouse104" instead"); it is no longer rewritten silently. Lookups still normalize their input.
- **Inputs:** `addSubject`/`addSession` refuse unknown fields (`unsupportedField`) instead of ignoring them.
- **UUID case:** every UUID comparison in the store is case-insensitive.
- `iCompensateFailedRebindDetach` is private.

### Format-neutral "is a dataset" (audit §12.1 item 5)
`subFunc/isUmitDataset(folder)`: a `.dat` with a valid UMIT header, or an `AcqInfos.mat` with a scalar `AcqInfoStream`. The store uses it for `addSession`, relocation and validation; `AcqInfos.mat` alone is no longer required (tested with a header-only folder). Errors: `Umitoolbox:UMITProjectStore:notADataset` (code `not_a_dataset`; replaces `missingAcqInfos`), `invalidAcqInfos` unchanged.

### Error identifiers added in Phase 6
`Umitoolbox:UmitFolder:invalidRootOverride`; `Umitoolbox:UMITProjectStore:invalidTarget`, `:unsupportedField`, `:notADataset`, `:bindingNotStale`, `:sessionAlreadyBound`, `:rawFolderMismatch`, `:rawFolderUnknown`, `:bindFolderFailed`. Removed: `:rigAssignmentNotSupported`, `:missingAcqInfos`, `:sessionRigMigrationFailed`.

## Phase 7 — Adoption (backend)

§4 Q7: a library match is shown only; otherwise the reference can be added with confirmation. `DataParams.registration` is never modified.
- `ProjectManagement/planReferenceAdoption(saveFolder)` (read-only): `notBound`, `notAligned`, `inLibrary` (the matching library resource, active or archived, by content hash) or `canAdd`, with a one-line message.
- `ProjectManagement/adoptFolderReference(saveFolder, 'Name', ..., 'Description', ..., 'SetActive', tf)`: only in state `canAdd`. It wraps the committed image in an ImageReference v2 with the bound project/subject/session identity (source type `adoptedFromFolder`), checks that the hash equals the commitment, and imports it with `addImageReference`. Otherwise `Umitoolbox:ProjectManagement:cannotAdoptReference`. A reference already in the library is never added twice.
- Tests `test/ProjectManagement/TestReferenceAdoption.m`: not bound / not aligned; add once (DataParams byte-identical, payload hash = commitment, second adoption refused); an existing (archived) library match is shown only, with nothing recorded (ledger, DataParams and library unchanged).
- The GUI offer is in DataViewer (`offerReferenceAdoption`, 2026-10-16); see "Pending GUI work".

## Phase 5–7 validation (2026-10-15/16, batch MATLAB through `validateTaskChanges`)

| Batch | Targets | Result |
|---|---|---|
| F | root isolation, `TestImportEmptyRigLibrary`, `TestUMITRigStore`, `TestUMITProjectStore` | 134/134 (the empty-library guards now run: 0 incomplete) |
| G | status helpers, `test/Analysis/ImageRegistration` (with the no-store scenarios), `test/ImageReference`, rewriters | 116/116 |
| H | `test/ProjectManagement`, `test/ImageReference`, `TestCreateRegistrationTform`, root isolation (after the store refactor) | 142/142 |
| I, J | binding states, stale links, transactional bind, RawFolder authority, adoption | 3 failures, all in new test fixtures (fixed) and one existing test whose flow item 13 now forbids (updated) |
| K1 | `test/RigManagement` (full rig suite), DataImport, ImagesClassification, importFromTif, HemoCompute, ImageRegistration, `test/subFunc`, DataParams, ImageReference | 356 run, 354 passed; 2 fixed and re-run in L |
| K2 | DataViewer, PipelineManager (both folders), EventsManager, Test_IO, Test_ROI, GUI theme, every other Analysis function folder | 700 run, 698 passed; 2 fixed and re-run in L |
| L | re-run of the fixed classes | 55/55 |

**Failures found by the regression runs:**
- `TestHemoCompute/testArchivedLibraryRigStillAnalyses` passed only because the real rig library held other rigs. On the isolated root its rig is the Active one, so archiving needs a replacement; the test now passes one when needed.
- `TestStreamingReadersHeaderedInput` (two HemoCompute cases) was already failing before this work. It bound its rig by stamping `AcqInfos.mat`, while `HemoCompute` has required a rig ledger since the rig redesign. It now binds through `UMITRigLedger`.

**Static analysis:** no new messages in the changed or new files. `UMITProjectStore.m` has 29 (30 before; one unused input fixed).

**Real stores untouched:** a listing of `Documents/LabeoTech/umIT` (paths, sizes, modification times) is identical before and after every batch above.
