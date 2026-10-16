# AUDIT — Project Management (read-only investigation)

**Audited:** branch `dev`, commit `b233fa61b04a9a1a41d99ec5a195fe22b936840d` ("refactor(analysis)!: absorb IOIAnalysis cores into native Analysis functions"). Two files were untracked when the audit started: `AUDIT_rig_management.md` and `RigManagement/UMITRigLedger.m`. During the audit, parallel rig work added `IMPLEMENTATION_NOTES_rig_ledger.md` and `test/RigManagement/TestUMITRigLedger.m`. In the notes, Phase 1 (ledger core) is done and Phases 2–6 are not started. Its "Project-side touchpoints (for Phase 6; not yet changed)" list (`IMPLEMENTATION_NOTES_rig_ledger.md:125-131`) matches R1–R5 and R9 in §11.1. No project code had been changed by that work at audit time.

**Date:** 2026-10-14.

**Scope:** project management (store, schema, binding, GUI tools, helpers) and every consumer of project data. Rig code is treated as moving ground: this audit lists where project code touches rigs but does not analyse rig behaviour (see `AUDIT_rig_management.md`).

**Conventions**
- **VERIFIED** means read in code (a file:line, or a method name for `.mlapp` files). **INFERRED** means deduced but not confirmed by reading every path or by running code. **UNKNOWN** means it could not be determined; these items are listed again in section 13.
- Nothing in the repository was modified apart from creating this file. No MATLAB code was executed and no tests were run.
- `.mlapp` files are zip archives. To read them, they were copied and unzipped into a session scratchpad outside the repo, and their code was read in targeted ranges. Cites such as `DataViewer.mlapp › ensureProjectBinding` refer to methods in those copies; the `.mlapp` binary has no stable line numbers.
- Not read, per `.claude/rules/mlapp-workflow.md`: all `*_exported.m` files (including `GUI/ProjectBindingManager_exported.m` and `GUI/DataViewer/ImageReferenceManager_exported.m`).
- Short forms used in citations:
  - `UPS` = `ProjectManagement/UMITProjectStore.m`
  - `SCH` = `ProjectManagement/getUMITProjectSchema.m`
  - `RIG` = `RigManagement/UMITRigStore.m`
- The reading was split across four parallel read-only passes:
  1. store internals;
  2. GUI and DataViewer;
  3. image alignment and references;
  4. pipeline, import and provenance.

  I re-read the claims that drive decisions myself: binding, validation, `addSession` preconditions, the alignment gates, `createRegistrationTform` resolution, the rig stamp on load, the PBM bind sequence, locks and `saveMatAtomic`.
- Sections 1–11 are descriptive only. Observations and suggestions are confined to section 12.

**Terminology.** The brief mentions a "session registry/ledger". The project code has **no separate ledger file or per-machine index**. "Registry" and "ledger" in `UPS` docstrings (e.g. `UPS:1412`, `UPS:1486`, `UPS:4133`) refer to the `subjectRegistry`, `sessionRegistry` and `resourceRegistry` arrays stored inside `project.mat` and `subject.mat` (section 3). VERIFIED. The only file called a "ledger" in the codebase is the rig-side `UMITRigLedger.json` (`RigManagement/UMITRigLedger.m:1-15`, untracked, in progress).

---

## 1. Overview

### 1.1 Modules and files

| Area | File | Role | Status |
|---|---|---|---|
| Store (sole writer) | `ProjectManagement/UMITProjectStore.m` (7019 lines, `classdef … < handle`) | Create/open/list projects; subjects, sessions; folder binding (`.umitlink`); managed resources (image references, registration transforms); renames, moves, removal; validation; locks | VERIFIED |
| Schema | `ProjectManagement/getUMITProjectSchema.m` (v1 only, `SCH:31-33, 41-352`) | Folder names, files, required fields, naming rules, resource types, editable fields, field descriptions | VERIFIED |
| Internal helpers | `ProjectManagement/_projectInternal/saveMatAtomic.m`, `computeFileChecksum.m` | Atomic MAT write; SHA-256 file checksum | VERIFIED |
| Root resolution | `subFunc/getUmitFolder.m` | `%USERPROFILE%\Documents\LabeoTech\umIT\<sub>` (`:40-46`) | VERIFIED |
| Project GUI | `GUI/ProjectManagerTool.mlapp` (PMT) | Browse/create/edit projects, subjects, sessions; bind raw folder; relocate/repair/remove/move/delete | VERIFIED |
| Binding GUI | `GUI/ProjectBindingManager.mlapp` (PBM) | Modal dialog from DataViewer: assign the loaded SaveFolder to project/subject/session, or rebind an orphan | VERIFIED |
| Session dialog | `GUI/promptAddSessionDialog.m` | Session-metadata form used by PBM. PMT has its own private copy (`ProjectManagerTool.mlapp › promptSessionInfo`) | VERIFIED |
| GUI helpers | `GUI/getSessionAcquisitionDatePrefill.m`, `GUI/resolveRigDisplayName.m`, `GUI/buildSessionRigChoices.m` | Date prefill from `AcqInfos.mat`; rig display names; `buildSessionRigChoices` has no production caller (test only) | VERIFIED |
| Main app | `GUI/DataViewer/DataViewer.mlapp` | Resolves the loaded SaveFolder's binding (`refreshProjectContextForCurrentData`), gates two tools (`ensureProjectBinding`), launches PMT/PBM, passes context to child tools | VERIFIED |
| Image references | `GUI/DataViewer/ImageReferenceManager.mlapp`, `genImageReferenceStruct.m`, `validateImageReferenceStruct.m`, `readImageReferenceSource.m` | Create, edit, select, activate and archive subject-owned image references in the store | VERIFIED |
| Manual alignment | `GUI/DataViewer/ImageAlignmentTool.mlapp`, `subFunc/applyImageAlignmentToFolder.m` | Interactive alignment to a selected reference; registers a session-owned registration transform | VERIFIED |
| Pipeline registration | `Analysis/ImageRegistration/createRegistrationTform.m`, `applyRegistrationTformOnFolder.m` | Estimate a transform against the subject's active reference (needs a project); apply a stored transform (no project) | VERIFIED |
| Folder provenance | `subFunc/{create,load,save,validate}DataParams.m` | `DataParams.registration` record (copy of tform, reference image, reference UUID/checksum) | VERIFIED |
| Tests | `test/ProjectManagement/TestUMITProjectStore.m`, `TestProjectManagerTool.m`, `test/ImageReference/TestImageReferenceStruct.m`, `test/Analysis/ImageRegistration/*`, `test/DataParams/TestDataParams.m` | Coverage of the above (not run) | VERIFIED (existence and test names) |
| Docs | `docs/documentation/userDocs/apps/project-session-integration.html`, `tools/project-management-folder-binding.html`, `tools/image-reference-manager.html`, `tools/manual_image_registration_tool.html` | User documentation of the project workflow | VERIFIED (existence; key statements in §11.6) |

### 1.2 Dependency map (VERIFIED unless marked)

```
getUmitFolder ──► UMITProjectStore ◄── getUMITProjectSchema
                     │   ▲   ▲
     saveMatAtomic ◄─┤   │   └── validateImageReferenceStruct (resource payload validator, SCH:169)
computeFileChecksum ◄┤   │
                     │   ├── ProjectManagerTool.mlapp ──► ImageReferenceManager.mlapp ('manage')
                     │   ├── ProjectBindingManager.mlapp ──► promptAddSessionDialog.m
                     │   ├── DataViewer.mlapp ──► PBM, PMT, ImageReferenceManager, ImageAlignmentTool,
                     │   │                         DataViewer_Coreg2Cams (store passed, unused)
                     │   ├── ImageAlignmentTool.mlapp ──► ImageReferenceManager ('select'), applyImageAlignmentToFolder
                     │   └── createRegistrationTform.m (static resolveProjectBinding)
                     └──► UMITRigStore.ensureDatasetRigAssociation (UPS:1235, UPS:1996)
```

- **No project dependency** in `classes/PipelineManager.m`, `classes/EventsManager.m`, `PipelineManagerTool.mlapp`, `DataViewer_EventsManager.mlapp`, `Analysis/DataImport/**`, or any `Analysis/**` function other than `createRegistrationTform`. VERIFIED by a grep for `UMITProjectStore`, `umitlink`, `projectUUID`, `sessionUUID`, `subjectUUID` and `ProjectBinding` across `classes/`, `Analysis/` and `subFunc/`, and in the extracted app code.
- Rig code has no dependency on project code (`AUDIT_rig_management.md` §9). VERIFIED there, not re-checked here.

---

## 2. Project data model

### 2.1 Entities

There are exactly three persisted entity types: **project**, **subject** and **session** (`SCH:114-129`). Two resource types hang off them (section 2.3), and a binding record lives in the data folders (section 4). VERIFIED.

There is **no group, cohort, condition or experiment entity** anywhere in the schema (`SCH:41-352`). VERIFIED.

**Project: `project.mat`, variable `ProjectInfo`** (required fields `SCH:190-204`)

| Field | Type / rule | Set by | Editable |
|---|---|---|---|
| `schemaVersion` | numeric, = 1 | `create` (`UPS:211`) | no |
| `projectUUID` | lowercase UUID (`java.util.UUID`, `UPS:6871-6875`) | `create` (`UPS:176, 212`) | no |
| `projectName` | non-empty text (`UPS:166-171`) | `create` | **no**: not in `editableFields.project` (`SCH:281-283`) |
| `description`, `institution`, `lab`, `principalInvestigator`, `experimenters`, `projectStartDate`, `projectEndDate`, `notes` | text, text list, or datetime (`iNormalizeMetadataValue`, `UPS:6616-6648`) | defaults (`UPS:6570-6603`) / input | yes (`SCH:281-283`) |
| `createdOn`, `modifiedOn` | datetime | store | no |
| `subjectRegistry` | struct array `{uuid,id,displayName,relativePath}` (`SCH:259-263`) | store | no |

**Subject: `subjects/<SubjectID>/subject.mat`, variable `SubjectInfo`** (required fields `SCH:206-230`)

- `schemaVersion`, `uuid`, `subjectID` (managed ID), `displayName` (defaults to the ID, `UPS:1045-1049`).
- Descriptive fields: `description`, `species`, `strain`, `genotype`, `sex`, `dateOfBirth`, `ageAtProjectEntry`, `ageReference`, `weight_g`, `housing`, `diet`, `healthStatus`, `earTag`, `rfid`, `notes`.
- `createdOn`, `modifiedOn`.
- `activeImageReferenceUUID` (pointer).
- `sessionRegistry` (same record shape as `subjectRegistry`).
- `resourceRegistry` (section 2.3).
- Editable: the descriptive fields and `displayName` (`SCH:284-288`). VERIFIED.

**Session: `subjects/<SubjectID>/sessions/<SessionID>/session.mat`, variable `SessionInfo`** (required fields `SCH:232-257`)

- Identity:
  - `schemaVersion`, `uuid`, `sessionID` (managed ID);
  - `subjectUUID` and `subjectID` (a parent copy);
  - `displayName`.
- Descriptive fields: `description`, `acquisitionDate` (datetime or NaT), `behavioralTask`, `stimulusNotes`, `pharmacology`, `surgery`, `virus`, `anesthesia`, `notes`.
- Folder locators:
  - `processedDataFolder` (SaveFolder, absolute canonical path, `UPS:1278, 6877-6892`) and `processedDataBindingUUID`;
  - `rawDataFolder` and `rawDataBindingUUID`.
- Rig mirror: `rigUUID`, `rigID`.
- `activeRegistrationTransformUUID`, `resourceRegistry`, `createdOn`, `modifiedOn`.
- Editable: the descriptive fields and `displayName`. Folder, binding and rig fields are explicitly excluded (`SCH:289-305`; enforced in `updateSessionMetadata`, `UPS:1906-1923`). VERIFIED.

**Required vs optional.**
- "Required" in the schema means *the field must exist*, not that it must be non-empty (`iValidateMetadataFields`, `UPS:6455-6472`). VERIFIED.
- Non-empty is enforced only for:
  - `projectName` (`UPS:169`);
  - managed IDs (`UPS:5170-5174`);
  - a session's SaveFolder (`session_save_folder_required`, `UPS:5713-5721`);
  - a session's rig pointer (`session_rig_required`, `UPS:5695-5699`).

  VERIFIED.
- Older metadata files missing newly optional fields are padded in memory on load (`iLoadMetadata` → `iAddMetadataDefaults`, `UPS:4872-4877`) and persisted on the next update. VERIFIED.

### 2.2 Identity and naming rules (VERIFIED)

- **UUIDs are identity.** Project, subject, session, binding and resource UUIDs are random v4 values from `java.util.UUID`, lowercased (`UPS:6871-6875`). Global UUID uniqueness is validated (`iValidateGlobalUUIDs`, `UPS:6312-6452`).
- **Managed IDs** (`subjectID`, `sessionID`) double as folder names. The rules are:
  - pattern `^[A-Za-z0-9][A-Za-z0-9_-]*$` (`SCH:135`);
  - at most 64 characters (`SCH:136`);
  - case-insensitive uniqueness (`SCH:137`; `iAssertUniqueRegistryID`, `UPS:5207-5231`);
  - Windows reserved names rejected (`SCH:141-146`) and reserved prefixes `.`, `__`, `.umit` rejected (`SCH:140`).
- Non-conforming IDs are **silently rewritten** with `matlab.lang.makeValidName` (`UPS:5176-5178`). Subject IDs are unique per project and session IDs per subject (`UPS:1055-1056`, `UPS:1232-1233`).
- **Project folder name** is built from the project name plus the full UUID (`iMakeProjectFolderName`, `UPS:6678-6699`; suffix `__<uuid>` assumed by `resolveProjectRoot`, `UPS:339`). Users cannot choose the project location (`UPS:100-101, 147-151`).
- **Display names** are free text and editable without renaming folders (`UPS:77`).

### 2.3 Resources

The schema defines exactly two resource types (`SCH:158-185`). VERIFIED.

| Type | Owner | Folder (`<State>` = `active`/`archive`) | File name | Payload contract |
|---|---|---|---|---|
| `imageReference` | subject | `subjects/<SubjectID>/references/image/<State>/` (`SCH:106-107, 165-167`) | `ImageReference_<yyyyMMddTHHmmss>_<uuid8>.mat` (`UPS:3089-3095`) | `.mat` with variable `ImageReference` passing `validateImageReferenceStruct` (`SCH:168-169`) |
| `registrationTransform` | session | `subjects/<SubjectID>/sessions/<SessionID>/transforms/registration/<State>/` (`SCH:110-112, 176-178`) | `SessionToReference_<ts>_<uuid8>.mat` | **none** (`requiredVariable=''`, `validatorFunction=''`, `SCH:179-180`) |

- **Resource record** (in the owner's `resourceRegistry`, `SCH:265-277`): `uuid`, `type`, `fileName`, `relativePath`, `displayName`, `description`, `createdOn`, `modifiedOn`, `status` (`available|active|archived`, `SCH:156`), `archivedOn`, `checksum` (SHA-256 of the file bytes, `computeFileChecksum.m:22, 46-57`; computed at `UPS:3112`), and `sourceFile` (the absolute path of the imported file). VERIFIED.
- Camera coregistration was moved out of the project to the rig store (`SCH:182-185`). VERIFIED.

### 2.4 Storage on disk (VERIFIED)

```
%USERPROFILE%\Documents\LabeoTech\umIT\projects\          (getUmitFolder.m:40-46; UPS:111)
  <projectName>__<projectUUID>\
    project.mat
    .umit\transactions\  .umit\recovery\  .umit\logs\project.log  .umit\project.lock\lockInfo.mat (transient)
    subjects\<SubjectID>\subject.mat
                       \references\image\{active,archive}\ImageReference_*.mat
                       \sessions\<SessionID>\session.mat
                                           \transforms\registration\{active,archive}\SessionToReference_*.mat
<any SaveFolder>\UMITProjectBinding.umitlink               (MAT format, variable ProjectBinding)
<any RawFolder>\UMITProjectBinding.umitlink                (only if bound and different from the SaveFolder)
```

- Folder lists come from `SCH:49-69, 93-100`. The log file name is `project.log` (`iAppendLog`, `UPS:6548-6566`).
- Strict "allowed children" rules make any unexpected file or folder inside managed locations a validation **error** (`UPS:2980-2984, 6267-6310`). OS incidental files are ignored (`SCH:150-151`).

### 2.5 Schema versions and migration (VERIFIED)

- Only schema **version 1** exists (`SCH:41-352`). Any other version raises `unsupportedVersion` (`SCH:349-351`).
- `ProjectInfo.schemaVersion` must equal the loaded schema version, else the error `schema_version_mismatch` (`UPS:2948-2955`).
- The binding struct has its own `version` (=1, `SCH:78`).
- **There is no migration code** between project schema versions. The only "migration" mechanisms are:
  - in-memory default padding of optional fields (`UPS:4872-4877`);
  - `migrateSessionRigAssignment` (`UPS:1975-2000`), which fills a blank session rig pointer. It has **no production caller**: grep over `.m` files and the extracted apps finds only the definition and tests.
  - `bindProcessedDataFolder` handling for "legacy schema-v1 projects" whose session has an empty `processedDataFolder` (`UPS:1362-1381`).

### 2.6 Fields that could later support group analysis (list only)

- Subject: `species`, `strain`, `genotype`, `sex`, `dateOfBirth`, `ageAtProjectEntry`, `ageReference`, `weight_g`, `housing`, `diet`, `healthStatus`, `earTag`, `rfid` (`SCH:212-224`).
- Session: `acquisitionDate`, `behavioralTask`, `stimulusNotes`, `pharmacology`, `surgery`, `virus`, `anesthesia` (`SCH:240-247`).
- **No field exists for experimental group, condition, treatment arm or cohort**, and there is no free-form key/value or tag field. VERIFIED (`SCH:190-257`).

---

## 3. Central storage, registry and the "ledger"

### 3.1 Location, format and contents (VERIFIED)

- **Location:** `getUmitFolder('projects')` → `%USERPROFILE%\Documents\LabeoTech\umIT\projects` (or `$HOME/...`) (`getUmitFolder.m:40-46`; `UPS:95-134`).
  - The folder is created if missing (`UPS:126-133`; `getUmitFolder.m:71-73`).
  - There is **no environment-variable or configuration override**: `getenv` is used only for USERPROFILE/HOME (`getUmitFolder.m:40-44`), and `iGetEnvironmentValue` (`UPS:6914`) is used only for the lock's user and host fields (`UPS:4919-4922`).
  - The help text says `<USER>/Documents/LabeoTech/Config/umIT` (`getUmitFolder.m:11`), which does not match the code.
- **Format:** MATLAB MAT files, one named struct variable per file (`SCH:64-75`), written with `saveMatAtomic` (`UPS:4771-4793`).
- **Registries:**
  - `ProjectInfo.subjectRegistry` and `SubjectInfo.sessionRegistry`: records `{uuid, id, displayName, relativePath}` (`SCH:259-263`; `UPS:6744-6752`). `displayName` is mirrored from the entity file on update (`UPS:1866-1867, 1937-1938`).
  - `resourceRegistry` inside each subject and session file (section 2.3).
- **There is no per-machine index, cache or session ledger:**
  - Project discovery is a `dir` scan that loads every `project.mat` (`resolveProjectRoot`, `UPS:315-393`; `listProjects`, `UPS:395-527`).
  - Folder→session lookup goes through the `.umitlink` in the folder (`resolveProjectBinding`, `UPS:676-743`; `findSessionByDataFolder`, `UPS:952-1019`) or through full scans of every `session.mat` (`iFindSessionPathMatches`, `UPS:4682-4736`; `iFindResourceLocations`, `UPS:3458-3510`).

### 3.2 Per-machine or shareable?

- **Per Windows user per machine.** VERIFIED (root derived from `%USERPROFILE%`, no override).
- **Data folders moved to another PC:**
  - The `.umitlink` travels with the folder (section 4). On a machine without that project, the binding resolves as **orphan** (`isOrphanProjectBinding`, `UPS:596-623`), and DataViewer reports `missingProject` (`DataViewer.mlapp › classifyProjectBindingFailure`). VERIFIED.
  - That machine may then delete the link (`removeOrphanProjectBinding`, `UPS:625-674`) or overwrite it (`ReplaceOrphanBinding=true`, `UPS:3753-3769`), both from the GUI (`DataViewer.mlapp › recoverMissingProjectBinding`; PBM rebind mode). VERIFIED.
  - INFERRED: when the folder returns to the original machine, its session there no longer matches the link.
- **Project folder copied to another PC:** it can be found and opened, since resolution is by UUID scan and tolerates renamed folders (`UPS:336-343`).
  - `SessionInfo.processedDataFolder` and `rawDataFolder` are **absolute paths** (`UPS:1278`), compared after case-folding on PC/Mac (`iNormalizeComparisonPath`, `UPS:6862-6869`).
  - A different drive letter or UNC mapping gives `bindingPathMismatch` (`UPS:719-726, 4669-4679`), so every folder must be relocated by hand (`relocateSaveFolder`). INFERRED from the comparison code.
- **Multiple users or machines on one project** (e.g. a project folder on a share):
  - Not supported by configuration, because the root cannot be pointed at a share. VERIFIED.
  - If it were reached anyway, the lock is a directory created atomically with `java.io.File.mkdir` (`UPS:4887-4912`), recording `operation`, `createdOn`, `processID`, `userName`, `hostName` (`UPS:4914-4924`).
  - There is no wait or retry. Stale-lock clearing is manual (`clearStaleLock`, `UPS:2835-2895`; age-based only, no PID or host liveness check) and has **no GUI or production caller** (grep). VERIFIED.
  - Readers take no lock. INFERRED from the getters, which do not call `iAcquireWriteLock`.

### 3.3 Missing, corrupt or inconsistent central records (VERIFIED)

- **`open`** (`UPS:255-313`):
  - errors if `project.mat` or its variable is missing, or the UUID does not match;
  - otherwise runs `validate('Mode','quick')` and sets `IsReadOnly = ~report.isValid` (`UPS:310-312`).
  - `validate` writes nothing to disk. The pass-1 grep of the validation helpers found no save, move, delete, mkdir or copy, and the top-level body (`UPS:2897-3020`) only updates `LastValidationReport`.
- **Every normal mutation** first calls `iAssertHealthyForMutation`, a **full** validation that includes SHA-256 of every resource file (`UPS:4955-4980`; checksums only in full mode, `UPS:6032-6046`). Any error blocks the mutation and flips the store to read-only. A second full validation runs after the write (`iAssertValidAfterMutation`, `UPS:5142-5151`).
- **Validation reads every bound data folder, in quick mode too** (`iValidateSessionNode`, `UPS:5701-5814`):
  - Missing external folder → **warning** `external_data_folder_missing` (`UPS:5753-5759`). The project stays writable.
  - Folder present but `AcqInfos.mat` missing or invalid → **error** (`UPS:5762-5777`).
  - Raw folder without `*.bin` or `*.tif(f)` → **error** (`UPS:5779-5786`).
  - Folder present but `.umitlink` missing → **error** `missing_project_binding` (`UPS:5789-5798`).
  - `.umitlink` unreadable or mismatched → **error** `invalid_project_binding` (`UPS:5801-5812`).
  - Session without a SaveFolder → **error** `session_save_folder_required` (`UPS:5713-5721`).
  - Session without a rig pointer → **error** `session_rig_required` (`UPS:5695-5699`). Only well-formedness is checked; the rig store is not consulted (`UPS:5680-5684`).
- **Consequence:** one session with a deleted `.umitlink`, a deleted `AcqInfos.mat`, a removed raw file or a blank rig field makes the **whole project read-only** at the next `open`. VERIFIED (`UPS:310-312`).
- **Narrow repair path:** `relocateSaveFolder`, `relocateRawDataFolder`, `repairSaveFolderBinding` and `removeSessionFromProject` bypass `iAssertWritable` and accept validation errors only if every error carries one of the codes `missing_project_binding`, `invalid_project_binding`, `incomplete_data_folder_binding`, `missing_acq_infos`, `invalid_acq_infos` or `missing_raw_data_files` **on the selected session** (`iAssertBindingMutationAllowed`, `UPS:5111-5140`). On success they reset `IsReadOnly=false` (e.g. `UPS:4205, 5096`). There is no comparable repair path for `session_rig_required` other than the uncalled `migrateSessionRigAssignment`.
- **Crash leftovers:**
  - Transactions are per-file atomic renames with compensating rollback (`saveMatAtomic.m:53-95`; e.g. `UPS:3812-3839`), not a journal. INFERRED: a crash mid-sequence can leave partial state.
  - Nothing replays or cleans up `.umit/transactions` or `.umit/recovery` on `open` (`UPS:255-313`).
  - Rename, move and rebind deliberately keep their recovery snapshot on failure (`UPS:2115-2124`).
  - Leftover recovery or transaction folders are inside `.umit`, an allowed child (`UPS:2980-2984`). INFERRED: they do not fail validation.

---

## 4. Folder binding (`UMITProjectBinding.umitlink`)

### 4.1 Contents and format (VERIFIED)

- File `UMITProjectBinding.umitlink` (`SCH:69`): a MAT file read with `load(...,'-mat')` (`UPS:577-578`), holding the variable `ProjectBinding` (`SCH:75`).
- Exact fields (`SCH:83-91`; written at `UPS:1290-1301, 3773-3784`):

  | Field | Content |
  |---|---|
  | `version` | 1 |
  | `bindingUUID` | new UUID per binding |
  | `projectUUID` | owning project |
  | `subjectUUID` | owning subject |
  | `sessionUUID` | owning session |
  | `folderRole` | `processedDataFolder` or `rawDataFolder` |
  | `createdOn`, `modifiedOn` | datetime |

- **No path, no human-readable ID or name, no rig information** is stored. `claude-context/architecture-reference.md:121` mentions "UUIDs + hints", but the writer adds no hints. VERIFIED.
- Reads are strict: unknown or missing fields, a wrong version, non-lowercase UUIDs or NaT dates are rejected (`iValidateProjectBindingStruct`, `UPS:6943-7008`).

### 4.2 Lifecycle of the file (VERIFIED)

| Operation | Effect on `.umitlink` | Cite |
|---|---|---|
| `addSession` | Creates the SaveFolder link (`folderRole='processedDataFolder'`). Refuses if any link file already exists, orphan or not (`bindingFileExists`), or if the path is already registered in this project | `UPS:1197-1212, 1290-1333` |
| `bindRawDataFolder` | If raw == SaveFolder, shares the existing link and writes no new file. Otherwise writes a second link in the raw folder (`folderRole='rawDataFolder'`) | `UPS:3729-3738, 3771-3825, 5035-5109` |
| `bindProcessedDataFolder` | Binds only if the session's `processedDataFolder` is empty (legacy); otherwise delegates to `relocateSaveFolder` | `UPS:1362-1381` |
| `relocateSaveFolder` / `relocateRawDataFolder` / `repairSaveFolderBinding` | Writes the link in the new folder. Deletes the old folder's link **only if the old folder is reachable** and the link belongs to this session, or is malformed | `UPS:3946-4229` (old link `:4118-4171`) |
| `unbindRawDataFolder` | Deletes the raw link (or rewrites a shared one); needs the folder to be available | `UPS:4467-4581` |
| `unbindProcessedDataFolder` | Always errors (`unbindSaveFolderDeprecated`): **a SaveFolder cannot be unbound** | `UPS:1407-1418` |
| `moveSessionToSubject` | Rewrites `subjectUUID` in each bound folder's link | `UPS:2399-2414` |
| `rebindSessionToProject` | Rewrites `projectUUID`/`subjectUUID` in a second phase; no production caller | `UPS:2454-2714` |
| `removeSessionFromProject` | Deletes the links of **reachable** folders; silently skips unavailable ones | `UPS:4231-4323` (skip `:4272-4274`) |
| `removeSubjectFromProject`, `deleteProject` | Preflight requires **every** bound folder to be online with a matching link, then deletes all links | `UPS:1498-1625, 4325-4435` |
| `removeOrphanProjectBinding` (static) | Deletes a structurally valid link whose `projectUUID` is not found under the local projects root; refuses if the project resolves | `UPS:625-674` |
| Write mechanics | `saveMatAtomic` to `.UMITProjectBinding.umitlink.<uuid>.tmp.mat`, then `movefile` onto the final name. `saveMatAtomic` itself uses another temp and `.bak` file in the same folder | `UPS:4601-4618`; `saveMatAtomic.m:53-95` |

`CLAUDE.md` names "ProjectBindingManager" as the owner of `.umitlink`. In code, `UMITProjectStore` is the only reader and writer (`UPS:84-85`). The PBM app only calls the store. VERIFIED.

### 4.3 Folder roles and what depends on them (VERIFIED)

- **Roles** (`SCH:79-81`):
  - `processedDataFolder`: the SaveFolder; required for every session.
  - `rawDataFolder`: optional; set only through PMT ("Bind/Change Raw Folder", `ProjectManagerTool.mlapp › onBindRawFolder`).
- **What depends on them:**
  - DataViewer accepts only `folderRole == 'processedDataFolder'`; any other role is classified `conflict` (`DataViewer.mlapp › refreshProjectContextForCurrentData`).
  - `createRegistrationTform` resolves the binding from the SaveFolder (`createRegistrationTform.m:423-430`).
  - Nothing outside the store and PMT reads `SessionInfo.rawDataFolder`. DataViewer, PipelineManager, EventsManager and Coreg2Cams use the folder-local `DataParams.folders.RawFolder` instead (`DataViewer.mlapp › resolveDataRawFolderForCurrentData`, `persistRawFolderToDataParams`). The two records are never reconciled: the store has no DataParams calls (grep).

### 4.4 Resolving a folder to project/subject/session (VERIFIED)

`UMITProjectStore.resolveProjectBinding(folder)` (`UPS:676-743`) works in this order:
1. Read and validate the link.
2. `open(projectUUID)`, which scans the projects root and quick-validates.
3. Load subject and session by UUID; the session's parent must be the bound subject.
4. `SessionInfo.<role>BindingUUID` must equal `bindingUUID`.
5. The folder's current absolute path must equal the registered path.

It returns a context with project, subject, session and rig values (`UPS:728-742`).

| Situation | Outcome | DataViewer state |
|---|---|---|
| No link | `readProjectBinding` error (`UPS:571-574`) | `unbound` |
| Malformed link | `bindingReadFailed` | `invalid`: alert only, no GUI repair or removal path |
| Project UUID not found locally | `resolveProjectFailed` "not found" (`UPS:380-384`) | `missingProject`: Retry / Reassign / Remove binding / Cancel (`recoverMissingProjectBinding`). Offered only from the Image Reference/Alignment gates |
| Project found, session or subject UUID mismatch | error | `conflict` |
| Folder moved or copied (path mismatch) | `bindingPathMismatch` (`UPS:719-726`) | `conflict`; fix via PMT › Rebind/Repair |
| Project invalid | resolves, store read-only | `readOnly` |

- DataViewer distinguishes `missingProject` from `conflict` by checking whether the error message contains "not found" (`DataViewer.mlapp › classifyProjectBindingFailure`). VERIFIED.
- The resolution chain runs on **every data load** and on every open of the Management menu or the Assign Rig submenu (`DataViewer.mlapp › loadDataSource`, `ManagementMenuSelected` → `refreshRigAssignmentMenu` → `refreshProjectContextForCurrentData`). VERIFIED.

### 4.5 Multiplicity, moves, unbinding (VERIFIED)

- **Several projects:** no. One link file per folder, and the store refuses to overwrite a non-orphan link (`UPS:1207-1212, 3753-3768, 4049-4060`). Within a project, a path may be registered to only one session (`iFindSessionPathMatches`, `UPS:1197-1203, 3740-3747`).
  - The cross-project check relies only on the link file. A project cannot see other projects' registered paths (INFERRED from `iFindSessionPathMatches` scanning only its own project).
- **Raw and processed in the same folder:** allowed; they share one link (`UPS:4999-5109`).
- **Move between projects:**
  - `rebindSessionToProject` exists (`UPS:2454-2714`) but has no production caller (grep). It is two-phase without cross-store atomicity (`UPS:2466-2496`) and gives the target a new subject UUID.
  - PMT offers only `moveSessionToSubject` within a project (`ProjectManagerTool.mlapp › onMoveSession`).
- **Unbind:**
  - A SaveFolder cannot be unbound while its session exists (`UPS:1407-1418`). The only ways to make a SaveFolder unbound are removing its session or subject, deleting the project (all of which delete links of reachable folders), or removing an orphan link.
  - A RawFolder can be unbound (`UPS:1399-1404`).

---

## 5. Project lifecycle

### 5.1 Operations (VERIFIED)

| Step | API | Key preconditions | GUI entry |
|---|---|---|---|
| Create project | `UMITProjectStore.create(struct('projectName',…))` (`UPS:136-253`) | Non-empty name; staged under `.__creating_<uuid>` then moved | PMT › `onNewProject`; PBM "Create new project" |
| Edit project | `updateProjectMetadata` (`UPS:1819-1849`) | Writable, healthy; `projectName` not editable | PMT metadata editor |
| Add subject | `addSubject` (`UPS:1021-1122`) | Writable, healthy; unique ID; unknown fields silently ignored (`UPS:6605-6614`) | PMT › `onAddSubject`; PBM "New subject" |
| Edit/rename subject | `updateSubjectMetadata`, `renameSubjectID` / `renameEntity` (`UPS:1851-1891, 2032-2130, 2216-2260`) | Rename moves the folder and rewrites relative paths in all registries | PMT editor |
| Add session | `addSession(subjectID, sessionInfo)` (`UPS:1124-1345`) | **Existing SaveFolder with `AcqInfos.mat` containing scalar `AcqInfoStream`** (`UPS:5274-5317`); folder outside the project root; no existing `.umitlink`; path not registered in this project; calls `UMITRigStore.ensureDatasetRigAssociation(saveFolder)` (`UPS:1235`), which **writes `AcqInfos.mat`** if it has no rig (`RIG:476-493`); rig fields cannot be passed (`UPS:1170-1177`) | PMT › `onAddSession` (`uigetdir` + `promptSessionInfo`); PBM "New session" (`promptAddSessionDialog`) |
| Edit/rename session | `updateSessionMetadata`, `renameSessionID` (`UPS:1893-1973, 2132-2214`) | Folder, binding and rig fields rejected | PMT editor |
| Bind raw folder | `bindRawDataFolder` / `relocateRawDataFolder` | Folder contains a top-level `.bin/.tif/.tiff` file (`UPS:5319-5333`) | PMT › `onBindRawFolder` |
| Relocate/repair SaveFolder | `relocateSaveFolder`, `repairSaveFolderBinding` | See §3.3 narrow repair path | PMT › `onRebindSaveFolder`, `onRepairSaveFolder` |
| Move session to subject | `moveSessionToSubject` (`UPS:2262-2452`) | Every bound folder must be reachable (`UPS:2369-2383`) | PMT › `onMoveSession` |
| Remove session | `removeSessionFromProject` (`UPS:1484-1496, 4231-4323`) | Session metadata and its resources are deleted permanently; external data is untouched | PMT › `onRemoveSessionFromProject` |
| Remove subject | `removeSubjectFromProject` (`UPS:1498-1566`) | All bound folders online with matching links | PMT › `onRemoveSubjectFromProject` |
| Delete project | `deleteProject` (`UPS:1568-1625`) | Same preflight; skips the writable and health checks; the Delete button is always enabled (`ProjectManagerTool.mlapp › configureActionButtons`) | PMT › `onDeleteProject` |
| Archive | **No archive state for projects, subjects or sessions.** Only resources have `available/active/archived` (`SCH:156`) | — | — |

### 5.2 What a user must do before processing data (VERIFIED)

- **Processing without a project:** import, view, pipeline (single and batch), events, ROI and mask, dual-camera coregistration, and `applyRegistrationTformOnFolder` need **no project step at all** (section 7).
- **To use Image Reference Manager, Image Alignment Tool or `createRegistrationTform`**, the user must:
  1. **Import** the data first, so that the SaveFolder contains a valid `AcqInfos.mat`. `addSession` refuses a folder without it (`UPS:1195, 5274-5317`). Sessions cannot be pre-registered before import. Conversely, `funcTemplateAcquisitionInitializer` refuses a SaveFolder that already contains a `.umitlink` (`Analysis/funcTemplateAcquisitionInitializer.m:149`). So the order is fixed: import first, then bind.
  2. Make sure a rig can be resolved for the dataset (`UPS:1235`; an error there aborts `addSession`).
  3. Create or select a **project** (PBM or PMT).
  4. Create or select a **subject**.
  5. Create a **session** pointing at the SaveFolder. This also writes the `.umitlink` (`UPS:1327-1333`) and may write the rig into `AcqInfos.mat`.
  6. For alignment: the project must be **valid** (writable), because `ensureProjectBinding('RequireWritable', true)` is used (`DataViewer.mlapp › ImageAlignButtonPushed`).
  7. For alignment or registration: the subject must have an image reference (step 8), and for `createRegistrationTform` it must be the **active** one.
  8. To create that image reference: open Image Reference Manager on a bound folder and create it from the current data or a file.
- **Shortcut:** in DataViewer, opening Image Reference Manager or Image Alignment Tool on an unbound folder launches PBM in `assign` mode, which bundles steps 3–5 in one dialog (`DataViewer.mlapp › ensureProjectBinding`, `launchProjectBindingManager`).
- **PMT launched from DataViewer** always opens the generic project chooser; the current binding is not passed in (`DataViewer.mlapp › openProjectManagerTool`).

---

## 6. Image alignment and image references (critical)

### 6.1 What an image reference is (VERIFIED)

- **Payload:** the struct `ImageReference`, version 1, with exactly 25 fields. Missing and unknown fields are both rejected (`validateImageReferenceStruct.m:74-87`; built by `genImageReferenceStruct.m:142-167`). The fields are:
  - `version`;
  - `image` (2-D finite `single`; non-finite values set to 0, `genImageReferenceStruct.m:139-140`);
  - `name`, `description`;
  - `projectUUID`, `projectName`, `subjectUUID`, `subjectID`, `sessionUUID`, `sessionID`;
  - `createdOn`, `lastModified`;
  - `imageSizeYX`, `axisConvention`;
  - `sourceFolder`, `sourceFile`, `sourceFrame`, `sourceType`, `processingNotes`, `createdBy`;
  - `wasManuallyTransformed`, `parentResourceUUID`, `manualTransformTform`, `manualTransformAppliedOn`, `manualTransformNote`.
- **Required non-empty:** `name`, `projectUUID`, `projectName`, `subjectUUID`, `subjectID`, `axisConvention`, `createdBy`. `projectUUID`/`subjectUUID` must match an RFC-4122 regex (`validateImageReferenceStruct.m:116-143, 221-244`).
- The payload deliberately carries **no** resource UUID, checksum or path (`genImageReferenceStruct.m:65-66`).
- **File:** a `.mat` file holding `ImageReference`, stored in the subject's `references/image/<State>/` folder in the central project (section 2.3). There is no separate image file.
- **Identity:**
  - The resource UUID is assigned at import (`UPS:3091`).
  - The SHA-256 checksum is of the **MAT file bytes** (`UPS:3112`; `computeFileChecksum.m`), not of the image content. INFERRED: re-saving an identical image gives a different checksum, because MAT headers embed a timestamp.
- **Link to the project:**
  - The registry record lives in the owner subject's `resourceRegistry`; the active pointer is `SubjectInfo.activeImageReferenceUUID` (`SCH:164, 228`).
  - The payload also embeds `projectUUID`/`subjectUUID`. The store does **not** check those against the owner at import (`iValidateManagedResourcePayload`, `UPS:3149-3212`, runs only the struct validator). Only `createRegistrationTform` checks them (`createRegistrationTform.m:496-501`).
- **Creation path** (`ImageReferenceManager.mlapp › createManagedImageReference`):
  1. `readImageReferenceSource` (TIF, `.dat` via `DatImageSource`, `.umt` via `UMTImageSource`, 2-D MAT variable; it averages up to 100 evenly spaced frames, `readImageReferenceSource.m:84-223`, with no project dependency).
  2. `genImageReferenceStruct` using the store's `ProjectInfo`/`SubjectInfo`.
  3. `save` to a temp file.
  4. `ProjectStore.addImageReference` (`UPS:1627-1636` → `iAddManagedResource`, `UPS:3033-3147`).
  - The first reference of a subject becomes active automatically; later ones are `available` (`UPS:3114-3121`).

### 6.2 Every project requirement on the alignment / reference path

| # | Location | What is required | Technically necessary or merely enforced |
|---|---|---|---|
| P1 | `DataViewer.mlapp › ImageAlignButtonPushed` | `ensureProjectBinding('RequireWritable', true)`: an unbound folder opens PBM; a `readOnly` project is refused | **Enforced.** The project is used only to pick the reference and to register the transform afterwards, which is optional (P4). VERIFIED |
| P2 | `DataViewer.mlapp › getImageAlignmentToolContext` / `getImageReferenceManagerContext` | Error `DataViewer:MissingProjectContext` unless project, subject and session UUIDs/IDs are all present (`hasProjectSessionContext`) | **Enforced.** VERIFIED |
| P3 | `ImageAlignmentTool.mlapp › SelectImageReferenceButtonPushed` | When `~HasProjectContext`, alert "No project context" and return; the tool cannot obtain any reference without a project | **Enforced.** Everything downstream uses only the image and display fields (P6). VERIFIED |
| P4 | `ImageAlignmentTool.mlapp › registerRegistrationTransformResource` | `addRegistrationTransform` + `setActiveRegistrationTransform`; returns `''` without context; failures downgraded to a warning | **Optional by design.** Its own help says the central store "must never block or roll back the local folder registration". VERIFIED |
| P5 | `DataViewer.mlapp › ImageRefButtonPushed`; `ImageReferenceManager.mlapp › startupFcn` | `ensureProjectBinding('RequireWritable', false)`; the manager errors `ImageReferenceManager:MissingProjectStore` / `MissingSubjectID` | **Necessary under the current storage model** (store-owned files). Not necessary to compute or hold the image. VERIFIED |
| P6 | `genImageReferenceStruct.m:107-114`; `validateImageReferenceStruct.m:116-143` | Non-empty, regex-valid `projectUUID`/`subjectUUID`, plus names | **Enforced.** No computation uses these fields except the identity check in P8. VERIFIED |
| P7 | `UPS:3036-3038` (`iAddManagedResource`) | Writable, healthy store and an owner subject | Necessary only for *managed* references. VERIFIED |
| P8 | `createRegistrationTform.m:418-503` (`iResolveManagedImageReference`) | `resolveProjectBinding` (link present, UUIDs consistent, path equal), then `getActiveImageReference(subjectID)`, then type/status/owner checks, then SHA-256 re-check, then payload validation, then payload UUIDs must equal the binding. The docstring says "Folder-local reference images are not used as a fallback" (`:47`) | **Enforced.** The registration math needs only the 2-D fixed image `refFr` (`createRegistrationTform.m:154-248`); there is no parameter to supply one. VERIFIED |
| P9 | `subFunc/applyImageAlignmentToFolder.m` | No project use; stores `opts.resourceUUID` only when given (`:1059-1064`); requires `AcqInfos.mat` (`:113, 313-331`) but does not use its content | Project: none. `AcqInfos.mat`: **enforced, not used.** VERIFIED |
| P10 | `Analysis/ImageRegistration/applyRegistrationTformOnFolder.m` | `DataParams.registration.tform` only | None. VERIFIED |
| P11 | `UPS:310-312` | An invalid project opens read-only, which blocks P1 (writable required) | Consequence of whole-project validation (§3.3). VERIFIED |

- Only two DataViewer callbacks call `ensureProjectBinding`: `ImageRefButtonPushed` and `ImageAlignButtonPushed` (grep of the extracted DataViewer code). VERIFIED.
- Both tool buttons stay enabled regardless of binding, so that their callbacks can launch PBM (`DataViewer.mlapp › updateGUIEnabledState`). VERIFIED.

### 6.3 What is recorded in the data folder after alignment

There are three write paths, and they produce different records. VERIFIED unless marked.

**(a) `createRegistrationTform`** (pipeline, estimate only; `createRegistrationTform.m:308-339`, saved through `validateDataParams` + atomic `saveDataParams`):
- `DataParams.registration` fields:
  - `isRegistered=false`, `isReviewed=false`;
  - `tform` (similarity), `transformType`, `method`;
  - `referenceDescription`;
  - `referenceFile` = **absolute path into the central project** (`:317`);
  - `referenceImage` = **full copy of the fixed image** (`:318`);
  - `imageReferenceUUID` and `imageReferenceChecksum` (the registry SHA-256) (`:319-320`);
  - `createdOn`, `source` (absolute moving-file path);
  - `qcStatus`, `qcWarning`, `qcMetrics` (MI before/after/delta, translation, rotation, scale, determinant);
  - `qcFigureFile`, `qcPreviewImageFile`;
  - `applied*` and `resourceUUID` cleared.
- Files: `registrationQC_<ts>.fig/.png` in the SaveFolder.

**(b) `applyRegistrationTformOnFolder`** (pipeline, destructive apply):
- Rewrites `.dat` and image `.umt` files in place.
- Sets `isRegistered`, `isReviewed`, `appliedOn`, `appliedBy`, `confirmationMode` (`:199-205`).
- No backup zip. Does not transform `.roi`, `DataParams.mask` or `DataParams.view` (no such code in the file).

**(c) `ImageAlignmentTool` → `applyImageAlignmentToFolder`** (GUI):
- Writes a backup zip `alignment_bkp_<ts>` in the SaveFolder (`:155-158`).
- Transforms `.dat`, image `.umt` and `.roi` files, and warps `DataParams.mask`.
- Rewrites `DataParams.view` size, origin, pixel size and axis convention (`:942-985`).
- Updates `registration` (`iUpdateRegistrationStruct`, `:1021-1065`):
  - `tform`, `method='ImageAlignmentTool'`;
  - `referenceFile` = **project-relative path**;
  - `referenceImage` (copy), `referenceDescription`, `source`, `applied*`, `manualAdjustment` (a field not in the DataParams schema);
  - `resourceUUID` only when non-empty.
- **`imageReferenceUUID` and `imageReferenceChecksum` are never written** (grep of `applyImageAlignmentToFolder.m` returns no match); the existing struct is copied forward (`:1025-1031`).
- Saved with a raw `save(...)` (`:226`), not `saveDataParams`.
- Writes `manualAlignmentReport_<ts>.mat` (`:230-233`).
- In the project: a `SessionToReference_*.mat` resource containing `tform`, `tformInfo`, `manualAdjustment`, `referenceFile` (relative), `movingSourceFile`, `createdOn`, **set active** before the folder apply runs (`ImageAlignmentTool.mlapp › ApplyRegistrationtoFolderButtonPushed` calls `registerRegistrationTransformResource` before `applyImageAlignmentToFolder`). That resource stores no reference UUID or checksum.

**dataHistory:**
- None of these functions writes `dataHistory.mat` itself (grep).
- When run as pipeline steps, PipelineManager records a standard entry (function name, version, parameters, inputs; `classes/PipelineManager.m:5733-5902`). INFERRED for these specific steps.
- The GUI path does not go through PipelineManager. VERIFIED.

**Is the folder record enough without the project?**
- **Re-applying or inspecting the transform:** yes. `tform` plus `referenceImage` or `view.imageSizeYX` is all `applyRegistrationTformOnFolder` needs (`:146, 361-384`). DataViewer's "Preview Dataset Registration Reference" reads `referenceImage` without a project (`DataViewer.mlapp › readExistingImageReference`). VERIFIED.
- **Verifying which reference was used:** no.
  - The checksum is of the central MAT file and cannot be recomputed from `referenceImage`. VERIFIED (§6.1).
  - The UUID is meaningful only inside the store.
  - The GUI path records no UUID or checksum at all, and may carry stale values from an earlier `createRegistrationTform` run. VERIFIED.
- **Reproducing the estimate after apply:** only through the GUI backup zip. The pipeline apply keeps no backup, and `source` is an absolute path. VERIFIED.

### 6.4 If the reference came from an arbitrary file (VERIFIED unless marked)

| Component | Effect |
|---|---|
| `readImageReferenceSource` | Already works on any supported file; no project. It may open `listdlg` for multi-entry inputs (`:411-466`) |
| `genImageReferenceStruct` / `validateImageReferenceStruct` | Cannot build or accept a struct without non-empty, valid `projectUUID`/`subjectUUID` and names |
| `UPS.iAddManagedResource` | Not usable without a store |
| `ImageAlignmentTool` | Blocked by P3. Past the guard it uses only `ImageReference.image`, `.name` and the reference file path fields |
| `applyImageAlignmentToFolder` | Does not validate the struct. Uses `referenceImage` / `.image` (`:92-97`), `.axisConvention` (`:971-975`), `.description` (`:1040-1044`). A plain image works |
| `createRegistrationTform` | Impossible today: no parameter for a reference image or file. Its `UseFile='auto'` mode also needs `ImageReference.sourceFile` to name a `.dat` present in the target folder (`:526-562`), so references built from TIF/MAT/UMT fail 'auto' even inside a project (`AutoSourceFileNotDat`) |
| `validateDataParams` | Already allows empty `imageReferenceUUID` and `imageReferenceChecksum` together (`:299-323`), so project-less provenance is schema-legal |
| `applyRegistrationTformOnFolder` | Unaffected |

Nothing outside `createRegistrationTform`, the DataParams helpers and tests reads `imageReferenceUUID` or `imageReferenceChecksum` (repo grep plus extracted apps).

### 6.5 Sharing, versioning, replacement, deletion

- **Across sessions of one subject:** yes; references are subject-owned (`SCH:161`). VERIFIED.
- **Across subjects:** no.
  - Image Reference Manager lists one subject and rejects other owners (`ImageReferenceManager.mlapp › selectReferenceRow`, `ResourceOwnerMismatch`).
  - `createRegistrationTform` rejects payload or owner mismatches (`:447-453, 496-501`).
  - Importing the same source for two subjects creates two independent copies with different UUIDs.

  VERIFIED.
- **Versioning:** no version numbers. "Edit" creates a new resource with `parentResourceUUID` set to the old one (`ImageReferenceManager.mlapp › applyReferenceEdit`). VERIFIED.
- **Replace:** `setActiveImageReference` moves the active pointer and demotes the previous resource to `available` (`UPS:1805-1810, 3233-3291`). VERIFIED.
- **Archive / restore:**
  - Archive moves the file to `archive/`, sets `archived`, and updates `relativePath`. Archiving the active resource requires a replacement UUID (`UPS:1649-1747`).
  - Restore returns it as `available`, not active (`UPS:1749-1803`).
  - The active pointer can never be cleared once set.

  VERIFIED.
- **Delete:** there is no resource-delete API. Resources disappear only with their owner (`removeSessionFromProject`, `removeSubjectFromProject`, `deleteProject`). VERIFIED.
- **Folders aligned with an archived or replaced reference:** there is **no staleness detection**.
  - No backlink exists from a reference to the folders aligned against it.
  - No code compares a folder's `imageReferenceUUID`/`Checksum` with the current reference.

  VERIFIED (grep). Consequences:
  - `createRegistrationTform` always re-reads the current active reference, so replacing it silently changes the next estimate. Asserted by the test `testNewActiveReferenceSupersedesFolderSnapshot`. VERIFIED.
  - The absolute `referenceFile` written by (a) points into `…/active/…` and dangles after archive. INFERRED from the paths.
  - A subject rename rewrites registry paths (`UPS:2081-2100`) but not paths already copied into DataParams. INFERRED.
  - Removing the subject deletes the reference file, while folders keep `referenceImage`. VERIFIED (removal deletes owner folders, `UPS:1527-1528`).
- **The ImageAlignmentTool selector** offers `available` references as well as the active one (`ImageReferenceManager.mlapp › refreshReferenceTable`, select mode). `createRegistrationTform` uses only the active one. VERIFIED.

### 6.6 Other project-owned resources (VERIFIED)

- **Registration transforms** (session-owned):
  - The only other project resource type. There is no payload contract (`SCH:179-180`).
  - Created only by `ImageAlignmentTool`.
  - **Nothing reads them back:** `getActiveRegistrationTransform` (`UPS:942-950`) has no app or `Analysis/**` caller (grep). Archive and restore work as for references, and there is no delete.
  - The folder copy of the transform (`DataParams.registration.tform`) is the one that is actually used.
- **Not project-owned:**
  - **Allen atlas:** `GUI/DataViewer/CCF_2022_mouse_topView_data.mat`, shipped beside the app (`AllenAtlasROICreator.mlapp › getDefaultAtlasInfoFilePath`).
  - **ROI files:** folder-local and schema-versioned (`subFunc/getROISchema.m`); no project fields.
  - **Masks:** `DataParams.mask`.
  - **`.umt` files.**
  - **Camera coregistration:** rig store (`SCH:182-185`).

---

## 7. Function requirements inventory

### 7.1 DataViewer-accessible tools and pipeline functions (VERIFIED unless marked)

| Function/Tool | File | Requires project? | Project data used | How it is obtained | Behaviour when no project |
|---|---|---|---|---|---|
| Data loading / viewing | `DataViewer.mlapp › loadDataSource`, `DatImageSource.m` | No | Binding context for the status label only | `refreshProjectContextForCurrentData` → `resolveProjectBinding` | Label shows "Not assigned"; everything works. (Side effect on load: rig stamp into `AcqInfos.mat`, §8) |
| Image Reference Manager | `ImageReferenceManager.mlapp` | **Yes** (read-only project allowed) | Store, ProjectInfo, SubjectInfo, session IDs, subject resources | Context from DataViewer (`getImageReferenceManagerContext`) or from PMT | Unbound: PBM opens in `assign` mode. `missingProject`: recovery dialog. `invalid`/`conflict`: alert, tool not opened. Launched without a store: startup error |
| Image Alignment Tool | `ImageAlignmentTool.mlapp` | **Yes (writable)** | Store, subject/session, the selected reference; registers a transform | `getImageAlignmentToolContext` | As above; `readOnly` refused with "Writable project required". The tool's own "no context" alert blocks reference selection |
| `createRegistrationTform` (pipeline) | `Analysis/ImageRegistration/createRegistrationTform.m` | **Yes** | Binding, active image reference, checksum, payload UUIDs | `UMITProjectStore.resolveProjectBinding(SaveFolder)` (`:423-430`) | Hard error `Umitoolbox:createRegistrationTform:BindingResolutionFailed` before any write (`:425-428`); test `testMissingBindingDoesNotMutateFolder` |
| `applyRegistrationTformOnFolder` (pipeline) | `Analysis/ImageRegistration/applyRegistrationTformOnFolder.m` | No | — | — | Works |
| `applyImageAlignmentToFolder` (backend) | `subFunc/applyImageAlignmentToFolder.m` | No (caller supplies the reference) | Optional `resourceUUID` | Argument | Works if given an image; requires `AcqInfos.mat` and `DataParams.mat` |
| Preview Dataset Registration Reference | `DataViewer.mlapp › readExistingImageReference` | No | — | `DataParams.registration.referenceImage` | Works |
| Dual-camera coregistration | `DataViewer_Coreg2Cams.mlapp` | No | Receives `ProjectStore`, subject/session IDs (stored but never used) and the **session** `rigUUID/rigID` | `OiSDUalCamCoregButtonPushed` passes them only if `hasProjectSessionContext` | Uses the default/active rig via `ensureRigContext` / `getOrCreateDefaultRig` instead of the session rig |
| Events Manager | `DataViewer_EventsManager.mlapp`, `classes/EventsManager.m` | No | — | — | Works |
| Pipeline Manager (single) | `PipelineManagerTool.mlapp`, `classes/PipelineManager.m` | No | — | SaveFolder/RawFolder from DataViewer (`DataParams.folders.RawFolder`) | Works |
| Batch pipeline | `DataViewer.mlapp › BatchPipelineMenuSelected`, `PipelineManagerTool.mlapp › selectDataFoldersInteractive` | No | — | User-built folder-pair list (`uigetdir`, parent-folder scan) | Works |
| Import functions | `Analysis/DataImport/ImagesClassification.m`, `importFromTif.m`, `Analysis/funcTemplateAcquisitionInitializer.m` | No | — | — | Work. The initializer refuses a folder containing `.umitlink` (`:149`) |
| HemoCompute | `Analysis/Imaging/HemoCompute.m` | No (needs a rig) | — | `AcqInfoStream.rigUUID/rigID` → rig store (`:370-405`) | Works; falls back to the default rig with a warning |
| ROI tools, Allen atlas, masks | `AllenAtlasROICreator.mlapp`, DataViewer ROI/mask methods | No | — | — | Work |
| Assign Rig menu | `DataViewer.mlapp › refreshRigAssignmentMenu`, `SelectRigMenuSelected` | No | Checkmark comes from the session `rigUUID` when bound | `ProjectContext` | Uses the dataset's `AcqInfos` rig |
| Project Manager | `ProjectManagerTool.mlapp` | Inherent | — | `listProjects` / `open` | Generic chooser; works without loaded data |
| Bind Current Session… | `DataViewer.mlapp › BindCurrentSessionMenuSelected` | Inherent | — | PBM `assign` mode | Requires loaded data |

### 7.2 Behaviour that differs when data is project-bound (VERIFIED unless marked)

- **Coreg2Cams:** project-bound mode passes the session-mirror rig; unbound mode falls back to the default/active rig, not the dataset's `AcqInfos` rig (`DataViewer.mlapp › OiSDUalCamCoregButtonPushed`; `DataViewer_Coreg2Cams.mlapp › ensureRigContext`).
- **Assign Rig checkmark:** session rig when bound, dataset rig when unbound (`refreshRigAssignmentMenu`). The context label shows the dataset rig first in both cases (`updateDataFolderContextLabel`).
- **Context label:** project, subject and session when bound, or a state message (`updateDataFolderContextLabel`).
- **`createRegistrationTform`:** works only when bound.
- **ImageAlignmentTool:** registers a session transform resource only when bound (P4).
- Nothing else in `Analysis/**`, `classes/**` or the pipeline changes behaviour with binding (grep).

---

## 8. Side effects and write-on-read

| # | Trigger (read-like action) | What is written | Where | Status |
|---|---|---|---|---|
| W1 | Opening any `.dat` in DataViewer (`DatImageSource` constructor, `DatImageSource.m:140-158`) | `AcqInfoStream.rigUUID`/`rigID` set to the **current Active Rig** when missing (`UMITRigStore.ensureDatasetRigAssociation` → `iSaveLoadedMatAtomic`, `RIG:476-493`) | SaveFolder `AcqInfos.mat` | VERIFIED |
| W2 | `addSession` (a project registration, not a read, but users may not expect it) | Same rig stamp, before the session is created; **not rolled back** if `addSession` fails later (`UPS:1235` vs `1246-1341`) | SaveFolder `AcqInfos.mat` | VERIFIED |
| W3 | Any project list, open or resolve (`getProjectsRoot`) | `mkdir` of the projects root if missing (`UPS:126-133`; `getUmitFolder.m:71-73`) | `%USERPROFILE%\Documents\LabeoTech\umIT\projects` | VERIFIED |
| W4 | DataViewer startup, file dialogs, Preferences dialog (`DataViewerPreferences.load`) | Creates `appPreferences.json` with defaults when missing (`GUI/private/readUmitPreferences.m:19-29`) | `getUmitFolder('preferences')` | VERIFIED |
| W5 | Every successful data open (`DataViewer.mlapp › rememberSuccessfulDataOpen` → `DataViewerPreferences.recordSuccessfulOpen`) | `appPreferences.json` re-saved unconditionally (`DataViewerPreferences.m:89-101`) | User preferences | VERIFIED |
| W6 | `PipelineManager.loadDataHistory`, including `isCachedResultOfPipeline`, whose help says it has no side effects (`classes/PipelineManager.m:4128` vs `4155`) | `dataHistory.mat` re-saved after pruning entries whose file is missing (`:5629-5655`; plain `save`) | SaveFolder | VERIFIED |
| W7 | Coreg2Cams launched on unbound data | `getOrCreateDefaultRig` may create a rig | Rig store | INFERRED from the name; rig behaviour out of scope |
| W8 | `DataViewer.mlapp › ensureDataParamsFileForCurrentFolder` | Creates `DataParams.mat` if missing | SaveFolder | VERIFIED that it creates the file; it is called from explicit actions (e.g. setting the RawFolder), not on plain load |

**Project-store reads that do not write (VERIFIED):**
- `open`, `validate` (quick and full), `listProjects`, `resolveProjectBinding`, `readProjectBinding`, `isOrphanProjectBinding`, `getSaveFolderBindingStatus`, `findSessionByDataFolder`.
- Optional fields are padded in memory only (`UPS:4872-4877`).
- No migration or stamping happens on open.

**Writes that look like reads but are user-confirmed:**
- `recoverMissingProjectBinding` → "Remove Binding" deletes the `.umitlink` after a confirmation dialog.
- PBM "Use Binding" writes nothing.

**Other project writes on mutation (expected):**
- `.umit/logs/project.log` append (`UPS:6548-6566`).
- Lock folder (`UPS:4880-4931`).
- Recovery and transaction folders (`UPS:5343-5359`).
- Temp files beside the `.umitlink` in data folders (`UPS:4610-4613`; `saveMatAtomic.m:53-54`).

---

## 9. Provenance and traceability

### 9.1 What data folders record (VERIFIED)

| File (owner) | Project-related content |
|---|---|
| `UMITProjectBinding.umitlink` (store) | `projectUUID`, `subjectUUID`, `sessionUUID`, `bindingUUID`, `folderRole`, dates. **No names, IDs or paths** |
| `DataParams.mat` (DataParams helpers) | No project, subject or session identity. `registration.imageReferenceUUID/Checksum` and `referenceFile` point into the central project; `registration.resourceUUID` points to a central session transform; `registration.tform` and `referenceImage` are full copies. `folders.RawFolder` (absolute) is a separate record from `SessionInfo.rawDataFolder`. `cameraCoregistration.*` points into the rig store (`createDataParams.m:114-238`; `validateDataParams.m:33-43, 222-244`) |
| `AcqInfos.mat` (importers, plus the rig store) | Raw acquisition facts (Width, Height, FrameRateHz, ExposureMsec, ImportedChannels, MultiCam, possibly `DateTime`) and `rigUUID/rigID`. **No project, subject or session identity.** `DateTime` is used only to prefill the session `acquisitionDate` (`GUI/getSessionAcquisitionDatePrefill.m:11-37`); whether `ReadInfoFile` produces it is UNKNOWN |
| `dataHistory.mat` (PipelineManager) | Per output file: function name, version, parameters, inputs (absolute SaveFolder/RawFolder paths), `pipelineVersion`. **No project, subject, session or rig identity** (`PipelineManager.m:5733-5902`) |
| `pipeLog.mat`, `events.mat` | No identity. `events.mat.RawFolder` is an absolute path (`EventsManager.m:1703-1719`) |
| `manualAlignmentReport_<ts>.mat`, `registrationQC_*.fig/.png`, `alignment_bkp_*` | Alignment details; the reference only as a relative path or image copy |

### 9.2 The deletion test: "if the project files were deleted, what could no longer be known or done?"

**No longer knowable from the data folders (VERIFIED from the schema; the "lost" conclusion is INFERRED by elimination over §9.1):**
1. **All project metadata:** name, description, institution, lab, PI, experimenters, start and end dates, notes (`SCH:190-204`).
2. **All subject metadata:** subject ID and display name (the link holds only a UUID), species, strain, genotype, sex, date of birth, age and age reference, weight, housing, diet, health status, ear tag, RFID, notes (`SCH:206-230`).
3. **All session metadata:** session ID and display name, description, acquisition date as entered (only `AcqInfos.DateTime` might remain), behavioural task, stimulus notes, pharmacology, surgery, virus, anesthesia, notes (`SCH:232-257`).
4. **Grouping:** which folders belong to the same subject or project. It is recoverable only by scanning folders and matching `subjectUUID`s in the links; the ID and name labels themselves are lost.
5. **The central RawFolder record.** The folder-local `DataParams.folders.RawFolder` survives.
6. **Image reference files** and their registry: status, active pointer, edit lineage (`parentResourceUUID`), archive history, display names and descriptions, original `sourceFile`.
7. **Registration-transform resources** (but the applied transform survives in `DataParams.registration.tform`).
8. **Verifiability of which reference a folder was aligned to.** `imageReferenceUUID`, `imageReferenceChecksum`, `referenceFile` and `resourceUUID` would point at nothing. Only the image copy `referenceImage` remains.
9. **The audit log** (`.umit/logs/project.log`).

**No longer doable:**
- `createRegistrationTform` (hard error).
- Image Reference Manager and Image Alignment Tool. Every link becomes an orphan, so they lead to the recovery dialog or PBM.
- Re-using a subject's reference for a new session.

**Still doable (VERIFIED):**
- Viewing, imports, single and batch pipelines including `applyRegistrationTformOnFolder`, EventsManager, ROIs, masks, HemoCompute, and dual-camera coregistration. HemoCompute and coregistration depend on the separate rig store.
- Previewing the stored registration reference.

### 9.3 Mirrors and duplicates between central records and folders (VERIFIED unless marked)

| Information | Central | Folder | Sync mechanism |
|---|---|---|---|
| Project/subject/session identity | `SessionInfo` / registries | `.umitlink` UUIDs | Written together by the store with rollback (§4.2); checked on every resolve and validate |
| SaveFolder path | `SessionInfo.processedDataFolder` (absolute) | — (the folder itself) | Path equality on resolve; manual relocate when moved |
| RawFolder path | `SessionInfo.rawDataFolder` | `DataParams.folders.RawFolder` | **None.** Two independent authorities (§4.3) |
| Rig | `SessionInfo.rigUUID/rigID` | `AcqInfoStream.rigUUID/rigID` | Copied once at `addSession` (`UPS:1235, 1281-1282`). `assignDatasetRig` updates only the folder. Validation never compares the two (`UPS:5680-5699`), so the mirror can go stale |
| Acquisition date | `SessionInfo.acquisitionDate` | `AcqInfoStream.DateTime` (if present) | Prefill only, at session creation |
| Reference image | Central `ImageReference_*.mat` | `DataParams.registration.referenceImage` (copy) | Copied at registration time; no later sync or staleness check (§6.5) |
| Registration transform | Session `SessionToReference_*.mat` (GUI path only) | `DataParams.registration.tform` | Written in the same GUI action; nothing reads the central copy |
| Image reference identity | Registry UUID + checksum | `DataParams.registration.imageReferenceUUID/Checksum` | Written by `createRegistrationTform` only; never by the GUI path |

---

## 10. Multi-session and batch features (VERIFIED unless marked)

- **Batch pipeline is folder-based, not project-based.** `PipelineManagerTool.mlapp › selectDataFoldersInteractive` builds SaveFolder/RawFolder pairs from `uigetdir` or a parent-folder scan with filters (`addPair`, `addSaveFolderOnly`, `addSaveFoldersFromParent`). `PipelineManager(saveFolders, rawFolders)` takes the list (`applyFolderPairTable`). There are no project, link or session references in this logic (grep).
- **Pre-filled first pair:** from the loaded file's SaveFolder plus `DataParams.folders.RawFolder` (`DataViewer.mlapp › BatchPipelineMenuSelected`), not from the project session.
- **Each folder runs independently.** There is no cross-session aggregation in `PipelineManager.executePipeline` (`classes/PipelineManager.m:1603-1673`), nor anywhere in `Analysis/**` / `subFunc/**` / `classes/**` (no subject, group or session identifiers used; grep).
- **Cross-session reuse of a reference:** sessions of the same subject share the subject's active image reference. `createRegistrationTform` can therefore run in a batch over sibling sessions, each resolving its own `.umitlink` (`createRegistrationTform.m:423-432`). This is the only multi-session feature that depends on the project.
- **Enumerating sessions by project:** only PMT's tree (`ProjectManagerTool.mlapp › refreshTree`) and the store's registries do this. No tool takes "all sessions of project/subject X" as input. There is no group concept (section 2.6).

---

## 11. Coupling with other subsystems

### 11.1 Rig management touchpoints (listed only; rig behaviour not analysed)

| # | Touchpoint | Location | Necessary or incidental |
|---|---|---|---|
| R1 | Session schema requires `rigUUID`/`rigID` | `SCH:252-253, 289-294` | Incidental (mirror of `AcqInfos`) |
| R2 | Validation error `session_rig_required` / `session_rig_incomplete` | `UPS:5680-5699` | Incidental policy; makes a no-rig session invalidate the whole project |
| R3 | `addSession` calls `UMITRigStore.ensureDatasetRigAssociation(saveFolder)`; rig fields in input rejected | `UPS:1170-1177, 1235, 1281-1282` | Incidental; writes `AcqInfos.mat` |
| R4 | `updateSessionMetadata` rejects rig fields | `UPS:1906-1923` | Incidental |
| R5 | `migrateSessionRigAssignment` → `ensureDatasetRigAssociation` + `iStampSessionRig` | `UPS:1975-2000, 4795-4838` | Incidental; no production caller |
| R6 | `resolveProjectBinding` returns `rigUUID/rigID` from `SessionInfo` | `UPS:736-737` | Incidental |
| R7 | DataViewer passes the session rig to Coreg2Cams when bound | `DataViewer.mlapp › OiSDUalCamCoregButtonPushed` | Incidental (the dataset holds the same information) |
| R8 | Assign Rig checkmark prefers the session rig when bound | `DataViewer.mlapp › refreshRigAssignmentMenu` | Incidental |
| R9 | PMT, PBM and `promptAddSessionDialog` show the dataset rig or Active Rig as a preview | `GUI/promptAddSessionDialog.m:14-27`; `ProjectManagerTool.mlapp › createRigDisplay`, `promptSessionInfo` | Display only |
| R10 | `GUI/buildSessionRigChoices.m` | `:44-46` | Test-only |
| R11 | New, untracked `UMITRigLedger.json` (folder-bound rig truth) and `_tasks/active/task_redesign_rig_management.md` ("Project Management out of scope, but must not break", per the pass-4 reading) | `RigManagement/UMITRigLedger.m:1-15` | Moving ground |

### 11.2 PipelineManager

- **No dependency on project code or structure.** VERIFIED (grep over `classes/`; `claude-context/status.md:131` agrees).
- The private `projectFolder` property (`PipelineManager.m:87, 167-173`) is a legacy name for the `.pipe` default location (`pwd`); it is unrelated to `UMITProjectStore`.
- The only project-dependent pipeline step is `createRegistrationTform`, which resolves the project itself from the SaveFolder.
- Coupling: none, so neither necessary nor incidental.

### 11.3 DataViewer

All dependencies listed below are VERIFIED.

| Coupling | Method(s) | Necessary? |
|---|---|---|
| Resolves the binding on every load and on Management/Assign Rig menu opens | `refreshProjectContextForCurrentData`, `classifyProjectBindingFailure`, `hasProjectSessionContext` | Incidental for display; necessary only for the two gated tools |
| Gates Image Reference Manager and Image Alignment Tool | `ensureProjectBinding`, `recoverMissingProjectBinding`, `launchProjectBindingManager` | Enforced (§6.2) |
| Builds child-tool context from `ProjectContext` | `getImageReferenceManagerContext`, `getImageAlignmentToolContext`, `OiSDUalCamCoregButtonPushed` | Enforced / incidental |
| Launches PMT and PBM | `openProjectManagerTool`, `BindCurrentSessionMenuSelected` | Inherent UI |
| Context label | `updateDataFolderContextLabel` | Display |

Startup takes no project arguments (`parseStartupArguments`, `startupFcn`).

### 11.4 EventsManager

- No project coupling (`classes/EventsManager.m`; `DataViewer_EventsManager.mlapp`; grep). VERIFIED.
- It reads `AcqInfos.mat` from the SaveFolder or the RawFolder fallback (`EventsManager.m:1795-1804, 2639-2670`).

### 11.5 Import functions

- No project code in `Analysis/DataImport/**` or `funcTemplateAcquisitionInitializer.m`. VERIFIED.
- Rig touchpoints (rig behaviour out of scope): `ImagesClassification.m:146, 206-208, 241-283`; `importFromTif.m:134-148, 337-339`.
- **Indirect project coupling:**
  - `funcTemplateAcquisitionInitializer` refuses a SaveFolder containing a `.umitlink` (`:145-153`). Incidental: a safety guard against overwriting a dataset.
  - Conversely, `addSession` requires the importer's `AcqInfos.mat` (`UPS:5274-5317`). This is necessary only because the store uses `AcqInfos.mat` as its "is a dataset" test.

### 11.6 Documentation statements (VERIFIED, paraphrased)

- `docs/.../apps/project-session-integration.html`: a project is optional for a single-recording pilot but highly recommended; a bound project is required by Image Reference Manager and Image Alignment Tool.
- `docs/.../tools/project-management-folder-binding.html`: a binding is required for the reference and alignment tools; users must not hand-edit `.umitlink`.
- `docs/.../tools/image-reference-manager.html`, `tools/manual_image_registration_tool.html`, `fcns/createRegistrationTform.html`: binding required.
- `claude-context/decisions.md:66-70`: a valid data folder must remain inspectable and processable without first registering it in a project.
- `claude-context/architecture-reference.md:121-128` describes the central store as authoritative for rig identity and mentions link "hints". Both are stale (§4.1, §11.1).

---

## 12. Observations and risks

*(Suggestions are allowed in this section only.)*

### 12.1 Relevant to making projects optional

1. **The alignment and reference path is the only real blocker, and most of it is enforced rather than necessary.**
   - The computations need a 2-D fixed image (P3, P8).
   - The folder already stores a full copy of that image plus the transform (§6.3).
   - `validateDataParams` already accepts empty reference UUID and checksum (§6.4).

   The gates to relax or bypass for a project-less mode are:
   - P1/P2: DataViewer `ensureProjectBinding` and context builders;
   - P3: `ImageAlignmentTool` selection guard;
   - P5: `ImageReferenceManager` startup;
   - P6: identity fields in the ImageReference struct;
   - P8: `createRegistrationTform`, which has no folder-local or file input.

   *Suggestion:* define a folder-local or file-based reference source (for example a reference struct whose project/subject identity is optional, or a schema v2) and make the store one of several sources rather than the only one.
2. **Checksums identify files, not images.** The current SHA-256 is over the central MAT file (§6.1) and cannot be checked against the folder's `referenceImage`. *Suggestion:* if provenance must end at the folder boundary, record a content hash of the image array itself (computed the same way for project and non-project references), so the folder can verify itself without the store.
3. **The GUI alignment path writes the weakest provenance.** It never writes `imageReferenceUUID/Checksum`, keeps stale values from an earlier estimate, uses a raw `save` that bypasses validation and atomicity, and adds an off-schema field (`applyImageAlignmentToFolder.m:226, 1021-1065`). The selected UUID and checksum returned by Image Reference Manager are dropped (`ImageAlignmentTool.mlapp › SelectImageReferenceButtonPushed`). Fixing this matters for any project-optional design, because the folder record becomes the only record.
4. **`.umitlink` is UUID-only.** A folder cannot say which subject or session it is without the store (§9.2). *Suggestion:* decide whether folder-local provenance should carry human-readable snapshots (subject ID, session ID, project name, perhaps selected subject attributes). This also serves the future high-throughput linking need.
5. **Fixed order: import before binding.** `addSession` requires `AcqInfos.mat`, and the initializer refuses folders containing a link (§5.2, §11.5). That suits optional projects (binding is a post-import step), but it ties the store to `AcqInfos.mat`, which the `.dat` header shift is changing. *Suggestion:* replace `iAssertValidSaveFolderDataset` with a format-neutral "is a dataset" test.
6. **The rig requirement leaks into project validity.** A session without a rig invalidates the whole project (`UPS:5695-5699`), and `addSession` writes `AcqInfos.mat` through the rig store (`UPS:1235`). With rig provenance moving to the folder ledger, the session rig mirror (R1–R8) looks like a candidate for removal rather than synchronisation. This is to be decided with the rig redesign.
7. **Two RawFolder authorities** (`SessionInfo.rawDataFolder` vs `DataParams.folders.RawFolder`) are never reconciled (§4.3). Only the folder-local one is used by processing. With optional projects, the central copy appears redundant.

### 12.2 Latent bugs and inconsistencies (code paths VERIFIED; outcomes INFERRED unless noted)

8. **Whole-project fragility.** Any one bad dataset (missing link, `AcqInfos.mat` or raw file, or a blank rig) makes the project read-only at `open` and blocks every normal mutation (§3.3). Every mutation also runs two full validations, including SHA-256 of all resources and I/O to every bound folder (`UPS:4966, 5145`). INFERRED: this scales poorly and fails on offline network drives only when the folder *exists but is incomplete* (a missing folder is only a warning).
9. **Stale links that no API can clear.** These cases leave a valid link to an *existing* project whose session no longer matches:
   - `removeSessionFromProject` skips offline folders (`UPS:4272-4274`);
   - relocation leaves the link in an offline old folder (`UPS:4118-4119`);
   - rebind phase 1 (`UPS:2563-2603`).

   `isOrphanProjectBinding` returns false for such a link (`UPS:610-613`), so `removeOrphanProjectBinding` refuses it and `addSession` fails with `bindingFileExists`. INFERRED: the folder then needs manual deletion of the link, which the docs forbid.
10. **Multi-machine hazard.** A link written by another machine's project looks like an orphan locally, and the GUI offers "Remove Binding" or "Reassign" (§3.2), which silently breaks the other machine's project. INFERRED.
11. **PBM bind is not transactional across entities.** `resolveBindingSelection` runs `create` → `addSubject` → `addSession` → `bindProcessedDataFolder` with no compensation in the caller's `catch` (`ProjectBindingManager.mlapp › resolveBindingSelection`, `BindFolderButtonPushed`). A failure leaves an empty project or subject. In "Create new project" mode a retry creates another project of the same name, because the new store is not kept. VERIFIED code path.
12. **PBM orphan rebind with a new session always fails.** The orphan link is still in the folder, and `addSession` refuses any existing link (`UPS:1205-1212`). The failure happens after the project or subject has already been created (item 11). VERIFIED code path.
13. **PBM rebind can take over a healthy session.** In `rebind` state, `resolveBindingSelection` skips the `sessionAlreadyBound` check, so a session correctly bound elsewhere can be relocated onto this folder. VERIFIED code path; INFERRED outcome (the original folder then resolves as `conflict`).
14. **DataViewer recovery gaps.**
    - "Bind Current Session…" always uses `assign` mode and cannot recover orphans.
    - The tooltip claims recovery options were shown at load, but `loadDataSource` never calls `recoverMissingProjectBinding`.
    - `invalid` and `conflict` states have no GUI repair path.
    - `missingProject` vs `conflict` is decided by message text ("not found").
    - No GUI exposes `clearStaleLock`.

    All VERIFIED in `DataViewer.mlapp › BindCurrentSessionMenuSelected`, `updateDataFolderContextLabel`, `ensureProjectBinding`, `classifyProjectBindingFailure`.
15. **PMT resource management is unreachable.** `refreshTree` never adds resource nodes, so Set Active, Archive, Restore and Open Location for resources (`onSetActiveResource` and related handlers) are dead. `onRetrySaveFolder` has no button. VERIFIED (`ProjectManagerTool.mlapp › refreshTree`, `configureActionButtons`). Image references are managed only through Image Reference Manager.
16. **Orphan active transform on a failed apply.** `ImageAlignmentTool` registers and activates the session transform *before* `applyImageAlignmentToFolder`. If the apply fails and rolls back, the project keeps an active transform that was never applied. VERIFIED order.
17. **Re-estimation erases the apply record.** `createRegistrationTform` unconditionally resets `isRegistered=false` and clears `applied*` and `resourceUUID` (`createRegistrationTform.m:308-339`). That re-enables `applyRegistrationTformOnFolder` on already-transformed data without `AllowReapply`. VERIFIED code; INFERRED double-transform risk.
18. **Two apply paths with different semantics.** The pipeline apply leaves ROIs, mask and view calibration untransformed and makes no backup; the GUI apply transforms them and keeps a backup (§6.3). VERIFIED.
19. **`referenceFile` has three meanings:** absolute central path, project-relative path, or source-file path in the dead `DataViewer.mlapp › startSetImageReference`. All of them can dangle after archive, rename or rebind (§6.5).
20. **`UseFile='auto'` (the pipeline default) fails for references not built from a `.dat`** (`createRegistrationTform.m:541-547`). VERIFIED.
21. **The store does not verify payload identity at import.** An image reference whose embedded project or subject UUID belongs elsewhere is accepted (`UPS:3149-3212`); only `createRegistrationTform` catches it. VERIFIED.
22. **`addSession` side effect is not rolled back:** the `AcqInfos.mat` rig stamp (`UPS:1235`) survives a later failure. VERIFIED.
23. **`deleteProject` cannot delete what its docstring promises.** It is documented for corrupted projects (`UPS:1577-1579`), but the preflight fails on any offline bound folder or unreadable session file (`UPS:4346-4411`). INFERRED: a crashed delete leaves `.deleteProjectStage_<uuid>`, which `resolveProjectRoot`/`listProjects` do not filter (`UPS:330-331, 424`).
24. **Smaller inconsistencies (VERIFIED):**
    - **IDs and inputs changed silently:** `makeValidName` rewrites IDs silently (`UPS:5176-5178`) and `addSubject`/`addSession` ignore unknown fields (`UPS:6605-6614`), while the update APIs reject them.
    - **UUID comparison case is inconsistent:** case-sensitive in `iAssertBindingMatchesSession` (`UPS:4627`) and resource lookups (`UPS:3475, 6545`); case-insensitive elsewhere (`UPS:796, 3622`).
    - **Public helper:** `iCompensateFailedRebindDetach` is public (`UPS:2716`).
    - **Store/test hygiene:** the `getUmitFolder` help path is wrong (`:11` vs `:46`); there is no root override, so tests write into the developer's real project store (`test/ProjectManagement/TestUMITProjectStore.m:139-142`; `test/Analysis/ImageRegistration/setupSyntheticRegistrationFolder.m:71`); see DFR-20261005-001 for the rig analogue.
    - **Duplicated and dead GUI code:** `promptAddSessionDialog.m` claims PMT as a caller but PMT has a private copy; `resolveRigDisplayName` is duplicated in DataViewer; `buildSessionRigChoices.m` is dead.
    - **Error-ID convention:** `readImageReferenceSource:*`, `ImageAlignmentTool:*` (inside a `subFunc`), `ImageReferenceManager:*`, `DataViewer:*` and `getUmitFolder:*` do not follow `Umitoolbox:<Area>:<condition>`.
25. **Robustness of writes.**
    - `saveMatAtomic` is a sound temp, verify, backup and rename sequence. Multi-file operations are compensating sequences, not journaled transactions, and nothing cleans up on `open`.
    - Readers take no lock.
    - The lock has no liveness check.
    - Absolute paths are not portable across machines or drive mappings. The same applies to `dataHistory.mat` input paths: moving a folder disables step skipping (`PipelineManager.m:5848-5852, 13995-14014`).

---

## 13. Open questions

1. **What does "project optional" mean for image references?** Options include: folder-local references only; a user-level reference library outside projects; or arbitrary files with content hashes. Who owns references created without a project, and can they later be adopted into a project?
2. **When a folder aligned without a project is later bound, should its existing `DataParams.registration` be accepted as is, re-verified, or linked to a project resource?**
3. **What should the folder carry for identity?** Should `.umitlink` (or a successor) carry human-readable subject and session IDs and selected metadata snapshots, so the folder survives loss of the central store (§9.2)? Which fields?
4. **Is a SaveFolder rig pointer still required for a session to be valid** (R2, R3) once `UMITRigLedger.json` lands, or should the session rig mirror be dropped?
5. **Should one dataset's problem keep invalidating the whole project** (§3.3, item 8), or should validation severity be per session?
6. **Is cross-machine or shared-project use in scope?** If yes, the root override, path portability and orphan-overwrite behaviour (item 10) need decisions.
7. **Which RawFolder record is authoritative:** `SessionInfo.rawDataFolder` or `DataParams.folders.RawFolder`?
8. **Should registration-transform resources stay?** They are write-only today (§6.6).
9. **Are `rebindSessionToProject`, `migrateSessionRigAssignment` and `clearStaleLock` intended to get GUI entry points, or are they candidates for removal?** None has a production caller.
10. **Which of the two apply paths (pipeline vs GUI) is canonical?** They differ for ROIs, mask, view and backup (item 18).
11. **UNKNOWN items not resolved in this audit:**
    - whether `ReadInfoFile` yields `AcqInfoStream.DateTime` for real `info.txt` files;
    - the exact `dataHistory` content for registration steps run through PipelineManager;
    - the runtime behaviour of `waitfor` on AppBase objects in `DataViewer.mlapp › openProjectManagerTool`;
    - whether co-located Raw==Save bindings can produce a link with `folderRole='rawDataFolder'`, which DataViewer would classify as `conflict`;
    - whether directory `movefile` is atomic on network shares or OneDrive (filesystem-dependent).
