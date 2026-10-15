# AUDIT — Rig Management (read-only investigation)

Date: 2026-10-14 · Branch: `dev` · Scope: rig management and every consumer of rig data. Project Management is out of scope except where flagged in section 9.

**Conventions**
- **VERIFIED** = read in code (file:line, or function name).
- **INFERRED** = deduced, not confirmed by reading every path or by running code.
- **UNKNOWN** = could not be determined; listed again in section 11.
- Nothing in the repository was modified. No MATLAB code was executed. `.mlapp` files are zip archives; to read them without touching the repo, I copied them to a scratchpad folder, unzipped them, and ran narrow greps on `matlab/document.xml`. Cites like `DataViewer.mlapp › refreshProjectContextForCurrentData` refer to methods in those extracted copies (no stable line numbers exist in the `.mlapp` binary).
- Not read, per `.claude/rules/mlapp-workflow.md`: all `*_exported.m` files.
- Not read: the `test/` suites (except grepping test names), `EventsManager` internals, and `docs/`.

**Terminology warning (read first).** The brief mentions a "centralized ledger" and a "binding file in save folders". Neither exists for rigs in this codebase:
- There is **no rig ledger**. The rig store is a folder tree plus one default-pointer file (section 3).
- There is **no rig binding file**. A save folder's rig association is two fields, `rigUUID` and `rigID`, stored inside `AcqInfos.mat` (section 4).
- The only "binding file" and "ledger" in the code (`UMITProjectBinding.umitlink`, the project session registry) belong to Project Management and carry no rig fields (section 4.6).

---

## 1. Overview

### 1.1 Modules and files

| Area | File | Role |
|---|---|---|
| Rig store (sole writer) | `RigManagement/UMITRigStore.m` (3687 lines, `classdef … < handle`) | Create/open/list/archive/activate rigs; spectra and filter-set library; coregistration resources; dataset↔rig association; optical resolution (`resolveOpticalConfiguration`) |
| Rig schema | `RigManagement/getUMITRigSchema.m` (v1, v2, v3; default v3 at L29-31) | Field lists, folder names, naming rules, spectrum grid |
| Built-in template | `RigManagement/resources/defaultRigs/OiS200.json` | The rig auto-created when none exists |
| Built-in library | `RigManagement/resources/library/{spectra/{camera,filter,illumination}/*.txt, filterSets/*.json, legacyMappings.json}` | Shipped spectra (16 files) and filter sets (`GCaMP`, `jrGECO`, `none`) |
| Legacy migration | `RigManagement/migrateLegacyOpticalLibrary.m`, `RigManagement/resources/legacyOpticalParams/*.mat` | One-time conversion of the old IOI MAT optics (not used at runtime) |
| Rig UI | `GUI/RigManagerTool.mlapp` | Create / duplicate / edit hardware / archive / restore / activate / validate rigs; coregistration files; spectra and filter-set library |
| Dataset association | `UMITRigStore.ensureDatasetRigAssociation`, `assignDatasetRig`, `readDatasetRigAssociation` (L446-575) | Reads/writes `rigUUID`/`rigID` in the save folder's `AcqInfos.mat` |
| DataViewer | `GUI/DataViewer/DataViewer.mlapp` (Rig menu, status label), `GUI/DataViewer/DatImageSource.m` (L140-158), `GUI/DataViewer/DataViewer_Coreg2Cams.mlapp` | Display, assign, load-time stamping; dual-camera coregistration |
| GUI helpers | `GUI/resolveRigDisplayName.m`, `GUI/buildSessionRigChoices.m`, `GUI/promptAddSessionDialog.m` | Display and preview only |
| Analysis consumers | `Analysis/Imaging/HemoCompute.m`, `Analysis/Imaging/+hemooptics/ioi_epsilon_pathlength.m`, `Analysis/DataImport/ImagesClassification.m`, `Analysis/DataImport/importFromTif.m` | The only `Analysis/**` functions that touch rig data (grep over `Analysis/`, `subFunc/`, `classes/`, `IOIAnalysis/` for `UMITRigStore`/`rigUUID`/`rigID`/`RigInfo`) |
| Local coreg provenance | `subFunc/createDataParams.m` L197-233, `loadDataParams.m` L123-170, `validateDataParams.m` L447-580 | `DataParams.cameraCoregistration` record (copy of tform + rig/resource IDs) |
| Provenance engine | `classes/PipelineManager.m` (`generateDataHistoryEntry` L5733, `isHistoryMatch` L13898) | `dataHistory.mat`, step skipping — contains no rig information |
| Project-side mirror | `ProjectManagement/UMITProjectStore.m`, `getUMITProjectSchema.m` | Stores a copy of `rigUUID`/`rigID` per session (section 9) |

### 1.2 Where the data lives (VERIFIED)
- Rig store root: `<USERPROFILE>\Documents\LabeoTech\umIT\rigs` (`subFunc/getUmitFolder.m` L40-46 builds `…\Documents\LabeoTech\umIT`; `UMITRigStore.getRigsRoot`, L92-124, appends `rigs` and `mkdir`s it if missing).
- Per-machine, per-Windows-user. No environment-variable override exists in `getUmitFolder`.

---

## 2. Rig model

### 2.1 Where rigs are defined
- A rig is a **folder** `rigs\<rigID>\` containing `rig.mat` (variable `RigInfo`, a scalar struct), plus subfolders (`UMITRigStore.create` L126-244; schema v3 at `getUMITRigSchema.m` L308-361).
- Folder layout (schema `folders`/`pathTemplates`, L183-212, v3 removes `calibration-files`):
  - `rig.mat`
  - `.umit/transactions/`, `.umit/recovery/`, `.umit/logs/operations.log`, `.umit/rig.lock/` (transient)
  - `transforms/camera-coregistration/{active,archive}/CameraCoregistration_<timestamp>_<uuid8>.mat` (file name built in `iAddManagedResource`, L2293-2297)
- Store-level (not per-rig) items under the rigs root:
  - `.umit/defaultRig.mat` — variable `DefaultRig` with `rigUUID` and `modifiedOn` (`iWriteDefaultRigUUID`, L3302-3313)
  - `.umit/store.lock/` — transient lock folder (`iAcquireStoreLock`, L3346-3365)
  - `library/spectra/<category>/<id>.txt` and `library/filterSets/<id>.json` — the user optical library (`iUserSpectrumFolder`, L3088-3092; `importFilterSet`, L902-928)

### 2.2 Fields of `RigInfo` (schema v3, `requiredRigFields`, `getUMITRigSchema.m` L339-353; `cameraFields`/`illuminationFields` L355-359)

| Field | Type | Used in computation? | Notes |
|---|---|---|---|
| `schemaVersion` | double (=3) | Validation only | `iValidateRigIdentity` L2837-2841 requires `== schema.version` |
| `uuid` | char, canonical UUID | **Identity** — the binding key | Immutable (`editableFields.rig` excludes it) |
| `rigID` | char `^[A-Za-z0-9][A-Za-z0-9_-]*$`, ≤64 | Identity/locator; is the **folder name** | **Renamable** (`renameRigID`, L1384-1455); also stored in each `AcqInfos.mat` |
| `displayName`, `description` | char | No (descriptive) | Editable |
| `createdOn`, `modifiedOn` | datetime | No | `modifiedOn` is the only change marker |
| `status` | `'active'`/`'available'`/`'archived'` | Lifecycle only | Exactly one `active` rig is required (`validateLifecycle`, L391-425) |
| `archivedOn` | datetime / NaT | No | |
| `metadata` | scalar struct (free-form) | **No consumer found** | Shallow-merged on update (`iApplyEditableUpdates`, L2678-2690). Only `create` and `duplicate` touch it. |
| `cameras` | struct array, ≤2 | **Yes** | Fields `index` (1 or 2, unique), `displayName`, `manufacturer`, `model`, `serialNumber`, `spectrumID` (`iValidateCameraRecords`, L2994-3031). Computation reads `index` and `spectrumID`. |
| `illuminations` | struct array | **Yes** | Fields `name` (`red`/`green`/`yellow`), `displayName`, `manufacturer`, `model`, `spectrumID` (`iValidateIlluminationRecords`, L3033-3071). Computation reads `name` and `spectrumID`. |
| `activeCoregistrationUUID` | char | **Yes** (at import and in Coreg2Cams) | Pointer to the active transform resource |
| `resourceRegistry` | struct array | Yes (resolves the coreg file) | Record fields (`resourceRecordFields`, L141-153): `uuid, type, fileName, relativePath, displayName, description, createdOn, modifiedOn, status, archivedOn, checksum, sourceFile` |

**Not in the rig.** The rig holds no frame rate, exposure, channel order, LED sequence/interleaving, binning, bit depth, pixel size, field of view, lens/magnification, filter set, or any acquisition timing (full-text check of schema v3 and `RigInfo` construction in `create`).

### 2.3 Spectra, filter sets, units (VERIFIED)
- **Spectra are not embedded in the rig.** A camera or illumination references them by `spectrumID` string. They live in two places, searched user-first then built-in (`iFindLibraryFile`, L3094-3102):
  1. User library: `<rigsRoot>\library\spectra\<category>\<id>.txt`
  2. Built-in library shipped with the code: `RigManagement\resources\library\spectra\<category>\<id>.txt`
- Categories: `illumination`, `filter`, `camera` (`schema.spectrum.categories`, L330).
- File format (e.g. `LED_521nm.txt`): a `#`-commented JSON metadata header (`id, legacyID, category, displayName, manufacturer, model, peakWavelengthNm, normalization`) followed by 301 rows of `wavelength<TAB>response`. Wavelength is in nm on a fixed 400:700 nm, 1-nm grid (`schema.spectrum.wavelengthNm`, L329). Response is unitless, normalized by its maximum into [0,1] (`getSpectrum` L674-683; `importSpectrum` L760-782 interpolates/normalizes on import).
- **Filter sets** (`library/filterSets/<id>.json`, built-in or user) have fields `id`, `displayName`, `excitationSpectrumID`, `emissionSpectrumID`. They are **not part of a rig**. They are chosen per analysis call via HemoCompute's `FilterSet` parameter. An empty ID or `'none'` resolves to a flat response of ones (`iResolveFilterSpectrum`, L3210-3214).
- Built-in filter sets: `GCaMP` (excitation = emission = `FF01496LP`), `jrGECO` (excitation empty, emission `FF01512630`), `none` (both empty). See observation O-14.
- Spectrum and filter-set **IDs are immutable by policy**: `importSpectrum` (L790-794) and `importFilterSet` (L918-921) refuse to replace an existing ID. `updateFilterSet` (L930-966) **can** overwrite a user filter set in place. Built-in files cannot be edited via the API, but they ship with the code and can change between releases.

### 2.4 Lifecycle: create, edit, store, archive, set active (VERIFIED)
- **Create**: `UMITRigStore.create(rigInfo)` (L126-244). Staged in a temp folder, atomically moved. The first rig created when no rig is active becomes `active`; otherwise `available` (L216-225). `getOrCreateDefaultRig` creates the built-in `OiS200` if the store is empty (L1088-1100).
- **Duplicate**: `duplicate` (L1125-1172) copies description/metadata/cameras/illuminations; gets a fresh UUID/timestamps; does **not** copy status, coreg registry or coreg files.
- **Edit in place** (all on non-archived rigs):
  - `setCameras` / `setIlluminations` / `setHardwareConfiguration` (L1174-1227) overwrite hardware definitions and bump `modifiedOn`.
  - `updateRigMetadata` (L1358-1382): displayName/description/metadata.
  - `renameRigID` (L1384-1455): renames the folder.
  - The RigManagerTool's "Apply Hardware Configuration" dialog offers "Apply correction" (in-place) or "Duplicate as new Rig" (`RigManagerTool.mlapp › onApplyHardware`).
- **Archive** (soft delete): `archiveRig(replacementRigUUID)` (L1527-1613). Archiving the active rig requires an `available` replacement, which becomes active in the same transition. `restoreRig` (L1615-1649) brings it back as `available` (or `active` if none is active).
- **Active**: `activateRig` (L1457-1525) demotes every other active rig to `available`; `setDefaultRig` (L1016-1030) is an alias for `activateRig` + writing the pointer. **"Default" and "Active" are the same thing** (the pointer is rewritten on every activation: L1477, L1522, L1609).
- **Immutability that exists today**: only *archived* rigs reject mutation (`iAssertHealthyForMutation`, L2619-2633, error `archivedRigReadOnly`). Active and available rigs are freely mutable.
- **Coregistration resources**: `addCameraCoregistration` / `importAndActivateCameraCoregistration` (L1651-1672) copy a `.mat` (must contain `tform`; optional `tformInfo`) into the rig with a checksum (`computeFileChecksum`, algorithm not inspected) recorded in the registry (`computeFileChecksum`; checked by `validate('full')`, L2142-2152). Resources can be archived/restored/purged (L1674-1876).
- **Versioning**: none. There is no rig version number, revision counter, or content hash for the rig record itself (only `modifiedOn` and a per-rig `operations.log`, section 7).

---

## 3. Ledger

**No centralized rig ledger exists.** (Searched for "ledger" across `*.m`, `*.md`, `*.json`: all hits are in `ProjectManagement/` and its tests — the *project* session ledger — plus a comment in `subFunc/appendImportedChannelInfo.m` that says it is *not* a ledger.)

What exists instead (all VERIFIED):

| Item | Location | Format | Written | Read |
|---|---|---|---|---|
| Rig records | `<rigsRoot>\<rigID>\rig.mat` | MAT, variable `RigInfo` | By every mutation (`iSaveRigInfo`, L2210-2217, via `saveMatAtomic`) | `open` (scans **every** rig folder and compares UUIDs: L255-287), `openByRigID`, `listRigs` |
| Default/active pointer | `<rigsRoot>\.umit\defaultRig.mat` | MAT, variable `DefaultRig = {rigUUID, modifiedOn}` | `activateRig`, `setDefaultRig`, `archiveRig` (with replacement), `getOrCreateDefaultRig` | `listRigs` (L384, **unguarded**), `getDefaultRig` |
| Per-rig audit log | `<rigRoot>\.umit\logs\operations.log` | TSV: `timestamp \t operation \t targetUUID \t result` (`iAppendLog`, L2706-2727) | Each mutation (errors silently swallowed) | Not read by any code found |
| Locks | `<rigRoot>\.umit\rig.lock\`, `<rigsRoot>\.umit\store.lock\` | Folder + `lockInfo.mat` | Acquire/release around mutations | — |

**How it is compared against bindings.** It isn't. A save folder's `rigUUID` is only looked up by `UMITRigStore.open(uuid)`, i.e. a directory scan of rig folders for a matching `RigInfo.uuid`. There is no stored list of "known rigs" to compare against.

**If the store/entry is missing** (VERIFIED by reading the callers):
- Rigs root missing: `getRigsRoot` silently creates an empty folder (L117-123). Then `open(uuid)` throws `Umitoolbox:UMITRigStore:openFailed` ("Rig UUID was not found", L289).
- Binding references a UUID not in the store: HemoCompute does not catch it, so the pipeline step errors (`HemoCompute.m` L370-371). `importFromTif` errors (L136-137, via `ensureDatasetRigAssociation` L467-472). DataViewer load **warns but continues** (`DatImageSource.m` L150-157 → `RigAssociationIssue`; `DataViewer.mlapp › (load path)` shows a "Rig Not Resolved" `uialert`).
- UUID unresolved but `rigID` present: **no fallback to `rigID`** in any of those paths — `rigID` is used only when `rigUUID` is empty (`HemoCompute.m` L370-373; `ensureDatasetRigAssociation` L467-481).
- `defaultRig.mat` malformed or pointing at a missing/archived rig: `listRigs()` throws (L384 reads the pointer without try/catch), which breaks `validateLifecycle`, `getActiveRig`, `listAssignableRigs`, the RigManagerTool list, and anything built on them. `getDefaultRig` throws `invalidDefaultRig` (L1006-1012).
- A single unreadable `rig.mat`: tolerated; shown as `IsReadable=false` in `listRigs` (L373-380) and skipped by `open`. `validateLifecycle` counts only readable rigs, so an unreadable active rig yields `zeroActiveRig`.
- Stale `store.lock`: `iAcquireStoreLock` fails if the folder exists (L3353-3356), with no age check and **no public API to clear the store-level lock**. `clearStaleLock(Force)` (L2011-2034) clears only a per-rig lock. Practical effect: every store-level operation, including `getOrCreateDefaultRig` (used at import time), fails until the folder is deleted by hand. (INFERRED consequence; the code path is VERIFIED.)

**Data moved to another PC without the rig store** (VERIFIED paths, section 5): the `AcqInfos.mat` carries `rigUUID`/`rigID`, but nothing carries the rig definition or spectra. The other machine's store has no matching UUID, so:
- DataViewer opens the data with a non-blocking warning;
- HemoCompute errors;
- importFromTif into that folder errors;
- Coreg2Cams ignores the dataset binding when not project-bound (section 5).

If that PC has its own rig with the same `rigID` but a different UUID, it is **not** matched (UUID takes precedence; `rigID` is only a fallback when UUID is absent).

---

## 4. Binding mechanism

### 4.1 What the "binding" is
- **Name/location**: no dedicated file. The binding is two fields, `rigUUID` and `rigID`, on the top-level struct inside `<SaveFolder>\AcqInfos.mat` (variable `AcqInfoStream`; `UMITRigStore.iLoadAcqInfoStream`, L2758-2788, also accepts `AcqInfos` or the first variable).
- **Format**: MAT-file; both fields are char vectors. No timestamp, no rig snapshot, no rig version, no author, no history.
- **Full rig-related content**: `AcqInfoStream.rigUUID`, `AcqInfoStream.rigID`. (`AcqInfoStream.ImportedChannels` carries `CamIdx`/`Color`/`DatFile`, which are *acquisition* metadata that rig resolution consumes but does not own.) A second, independent record exists in `DataParams.mat → cameraCoregistration.rigID/resourceUUID` (section 7).

### 4.2 Who creates it
| Code | When | Source of the rig |
|---|---|---|
| `ImagesClassification` (default `ApplyRig=true`) | After classification (L206-208) | `getOrCreateDefaultRig()` — **always the machine default**, even if the folder previously pointed to another rig (L146-160) |
| `importFromTif` | End of import (L337-339) | If the existing `AcqInfos.mat` already has a rig: that rig (error if archived, L138-141). Else `getOrCreateDefaultRig()` (L148) |
| `DatImageSource` constructor | Every time a `.dat` is opened in DataViewer, **if** `AcqInfos.mat` exists (L140-158) | `ensureDatasetRigAssociation` → the current **Active** rig (L480) when the folder has no binding |
| `UMITProjectStore.addSession` / `migrateSessionRigAssignment` | Project operations | `ensureDatasetRigAssociation` (`UMITProjectStore.m` L1235, L1996) |
| DataViewer › Rig › Assign Rig (`SelectRigMenuSelected`) | User action | `UMITRigStore.assignDatasetRig` (L500-544) |

### 4.3 Who reads it
`readDatasetRigAssociation` (L546-575, read-only, swallows all errors): DataViewer (`refreshProjectContextForCurrentData`, `refreshRigAssignmentMenu`, `SelectRigMenuSelected`), `ProjectManagerTool`, `promptAddSessionDialog` (L15-16). `HemoCompute` reads `AcqInfoStream.rigUUID`/`rigID` itself (L360-373). `ensureDatasetRigAssociation` reads and may write.

### 4.4 Can it change after creation? What happens to outputs?
- **Yes.** `assignDatasetRig` (L500-544) overwrites both fields with any active/available rig (archived targets rejected, L530-534). No confirmation is enforced in the backend (the DataViewer menu adds a `uiconfirm`).
- It rewrites only `AcqInfos.mat` (atomic replace via `iSaveLoadedMatAtomic`, L2790-2821). It does **not** touch existing outputs (`HbO.dat`/`HbR.dat`, coregistered `.dat` files, `dataHistory.mat`, `DataParams.cameraCoregistration`), does not log the previous rig, and does not mark anything stale.
- `ImagesClassification` and `importFromTif` also write/refresh it as a side effect of import.

### 4.5 Behaviour when the binding is missing
| Consumer | Behaviour |
|---|---|
| `DatImageSource` | Stamps the current Active rig into `AcqInfos.mat` (silent write at file-open time; `ensureDatasetRigAssociation` L488-494). Failure → warning + `RigAssociationIssue` (L150-157), load continues |
| `HemoCompute` | Uses the default/active rig with warning `Umitoolbox:HemoCompute:MissingRigUsingActiveRig` (L375-380); errors `MissingRig` if no single active default exists (L406-408). Does not write the binding |
| `importFromTif` / `ImagesClassification` | Binds the machine default (creating `OiS200` if the store is empty) |
| `promptAddSessionDialog` | Previews the Active rig (L17-26) |
| `readDatasetRigAssociation` | Returns empty strings |

### 4.6 The project binding (out of scope, noted for completeness)
`UMITProjectBinding.umitlink` (variable `ProjectBinding`) has fields `version, bindingUUID, projectUUID, subjectUUID, sessionUUID, folderRole, createdOn, modifiedOn` (`getUMITProjectSchema.m` L83-91). It contains **no rig fields**.

---

## 5. Rig resolution path (critical)

### 5.1 Trace, per consumer of rig data

| Consumer | Where the rig comes from | Uses folder binding? | Uses "active/default" rig? |
|---|---|---|---|
| `HemoCompute` (`localResolveOpticalInfo`, L350-384) | `AcqInfoStream.rigUUID` → `UMITRigStore.open`; else `rigID` → `openByRigID`; else `getDefaultRig`/single active (L386-409) | **Yes (primary)** | Only as fallback when the binding is empty, with a warning |
| `HemoCompute 'preflight'` (L332-348, called by PipelineManager) | Same `localResolveOpticalInfo` | Yes | Same fallback |
| `ImagesClassification` (L145-182) | `getOrCreateDefaultRig()` | **No** — ignores any existing folder binding | **Yes** (and may auto-create) |
| `importFromTif` (L133-155) | Existing binding if present, else `getOrCreateDefaultRig()` | Yes if present | Fallback |
| `DatImageSource` (L140-158) | Stamps the Active rig if no binding | Preserves existing | Yes, to fill a blank |
| `DataViewer_Coreg2Cams` (`ensureRigContext`, `openRigStore`) | `RigUUID`/`RigID` passed by DataViewer **only when project-bound** (`DataViewer.mlapp › OiSDUalCamCoregButtonPushed`: rig args built inside `if app.hasProjectSessionContext()`); else `getOrCreateDefaultRig()` | **No** — it never reads the dataset's `AcqInfos.mat` rig | **Yes** when not project-bound |
| DataViewer labels (`refreshProjectContextForCurrentData`, `updateDataFolderContextLabel`, `getCurrentRigStatusText`) | Dataset binding first, session mirror as fallback | Yes | No |

### 5.2 Direct answer: would analysing data bound to an archived rig work today?
**Yes for HemoCompute (VERIFIED by code, no dedicated test located).**
- `UMITRigStore.open(uuid)` does not filter by status (L246-290).
- `resolveOpticalConfiguration` documents "Archived Rigs remain valid for analysis of historical data" (L1232-1235) and never checks status.
- Searched `test/Analysis/HemoCompute/*` for "archiv": no match. `test/RigManagement/TestUMITRigOptics.m:291` covers only that an archived rig rejects *mutation*. This end-to-end path appears untested.

**Not everywhere:**
- `importFromTif` into a folder bound to an archived rig → error `ArchivedRig` (L138-141) (intended: archived rigs cannot receive new data).
- `ImagesClassification` ignores the folder binding and rebinds to the machine default.
- `assignDatasetRig` refuses an archived target (L530-534).
- `Coreg2Cams`: reading an archived rig's active coregistration can work, but **saving** a new coregistration to it fails (`archivedRigReadOnly`, L2629-2632, through `iAssertHealthyForMutation`). When not project-bound, the tool uses the machine default rig instead.

**Caveat that applies to all "bound" resolution.** The binding is a *reference* (UUID). Resolution reads the rig's **current** on-disk state (cameras, illumination spectrumIDs) and the **current** content of the referenced spectra files. If an active/available rig was edited in place after the data was acquired or analysed, the same binding silently resolves to different optics.

---

## 6. Function requirements inventory

Scope: every `Analysis/**` function (all are pipeline-compatible and reachable from DataViewer's pipeline tools) plus the DataViewer tools that touch rigs. Rig use was determined by grep for `UMITRigStore|rigUUID|rigID|RigInfo|resolveOpticalConfiguration` over `Analysis/`, `subFunc/`, `classes/`, `IOIAnalysis/`; `.mlapp` apps via extracted copies. Every function not listed with rig fields had **zero** matches.

### 6.1 Functions that touch rig data

| Function | File | Requires rig? | Rig fields used | How rig info is obtained | Behaviour when rig is missing |
|---|---|---|---|---|---|
| `HemoCompute` | `Analysis/Imaging/HemoCompute.m` | **Yes** | `cameras(index).spectrumID` (by each channel's `CamIdx`), `illuminations(name).spectrumID`, `uuid`, `rigID`; plus filter set `FilterSet` | Folder binding in `AcqInfos.mat`; fallback = default/active rig (L350-409). Calls `rigStore.resolveOpticalConfiguration` (L382) | Binding present but unresolvable → **error** (`openFailed`). No binding + no single active default → **error** `MissingRig`. No binding + active default → **silent-ish default** (warning only) |
| `ImagesClassification` | `Analysis/DataImport/ImagesClassification.m` | **Yes by default; optional** via `ApplyRig=false` (not exposed to PipelineManager, L61-67) | `uuid`, `rigID`; `activeCoregistrationUUID` → coreg `.mat` | `UMITRigStore.getOrCreateDefaultRig()` (L146) | None exists → **auto-creates** `OiS200` and warns (L154-160). Invalid active coreg resource → **error** (L174-181). No active coreg → imports without it (`noActiveResource`) |
| `importFromTif` | `Analysis/DataImport/importFromTif.m` | **Yes, no opt-out** | `uuid`, `rigID`, `status` | Existing binding, else `getOrCreateDefaultRig()` (L133-155) | None exists → **auto-creates** `OiS200` (no warning, unlike ImagesClassification). Bound to archived → **error** |
| `ioi_epsilon_pathlength` | `Analysis/Imaging/+hemooptics/` | Indirect (takes the `opticalInfo` struct) | Resolved spectra: `illuminationResponse`, `cameraResponse`, `excitationResponse`, `emissionResponse`, `activeRows` | Via HemoCompute | Asserts on missing/invalid fields (L35-53) |
| `applyTform2Cams` (`subFunc`) | `subFunc/applyTform2Cams.m` | No (receives `tform`) | — | Callers pass the rig's coreg payload | — |
| `DataViewer_Coreg2Cams` | `GUI/DataViewer/DataViewer_Coreg2Cams.mlapp` | **Yes** (hard) | `uuid`, `rigID`, `activeCoregistrationUUID`, `resourceRegistry` | Project-session mirror when project-bound; else `getOrCreateDefaultRig()` | Alert "No rig available"; cannot save |
| `DataViewer` (shell) | `GUI/DataViewer/DataViewer.mlapp` | **No** | `displayName`, `rigID`, `uuid` for labels/menu | `readDatasetRigAssociation`, `listAssignableRigs` | Label "Rig: not assigned"; `DatImageSource` warns but loads |

### 6.2 Functions with **no** rig dependency (zero matches)
`correctMotionArtifact`, `run_BloodFlow`, `ConvertToTiff`, `getEvents`, `GSR`, `apply_aggregate_function`, `normalizeBSLN`, `normalizeLPF`, `normalizeZScore`, `spatialGaussFilt`, `applyRegistrationTformOnFolder`, `createRegistrationTform`, `Ana_Speckle`, `HemoCorrection`, `SpeckleMapping`, `apply_detrend`, `genAmplitudeMaps`, `split_data_by_event`, `genRetinotopyMaps`, `genVSM`, `calculateResponseFeatures`, `genCorrelationMatrix`, `genImageTimeSeriesUMT`, `getDataFromROI`; the four `funcTemplate*` files; `DataViewer_EventsManager`, `ImageAlignmentTool`, `ImageReferenceManager`, `AllenAtlasROICreator`, `PipelineManagerTool`, `ProjectBindingManager` (extracted-`.mlapp` grep: 0 matches each). The `funcTemplate.m` header states the rule explicitly: frame rate and similar come from the `.dat` header or explicit arguments, "Never read the frame rate from AcqInfos.mat" (L39-43).

`createRegistrationTform` requires a **project** binding and project image reference (`UMITProjectStore.resolveProjectBinding`, L424), not a rig. Out of scope.

### 6.3 Hidden dependencies (implicit acquisition assumptions not declared in the rig)

1. **Channel-name/file convention.** `HemoCorrection` hard-codes reference channel file stems `red`, `green`, `yellow` (L119-129, parameter named `Amber`). `HemoCompute` maps illuminations to the fixed row order red/green/yellow (`ioi_epsilon_pathlength.m` L16-19, schema `canonicalIlluminations`). The rig defines which colors exist but not which `.dat` they correspond to; that mapping lives in `AcqInfoStream.ImportedChannels` (`Color`/`Tag`/`DatFile`/`CamIdx`), which `HemoCompute` consumes through `resolveOpticalConfiguration` (L1284-1355).
2. **Legacy inference.** If `ImportedChannels` is absent, `resolveImportedChannelFallback` infers channel mapping from `.dat` headers, legacy `IlluminationN` declarations, or **canonical filenames** (`subFunc/resolveImportedChannelFallback.m` header, L1-18; called by `HemoCompute.m` L367).
3. **Camera model ↔ rig.** `ImagesClassification` reads the camera model from `info.txt` and uses `D1024`/`D1312`/`BFLY` to choose dropped-frame logic (L787-808), defaulting to `'D1024'` with a warning when absent (L488-493). The rig's camera `model`/`spectrumID` is not consulted.
4. **CamIdx.** `HemoCompute` requires each channel's `CamIdx` ∈ {1,2} and a matching rig camera `index` (L1317-1332). Both import paths default a missing `CamIdx` to 1 (`importFromTif.m` L324-328; `ImagesClassification.m` L684-686).
5. **Frame rate/exposure.** `HemoCompute` reads per-channel `frameRateHz`/`exposureMsec` from `.dat` headers (L458-521), not from the rig; it enforces equal durations (1e-3 s tolerance) and refuses to upsample. `Ana_Speckle`, `run_BloodFlow` and others resolve exposure from the `.dat` header (`resolveDataInfoValue`) with explicit-override.
6. **Interleaving (LED sequence).** Captured at import into `AcqInfoStream.IlluminationN` and `ImportedChannels`; not in the rig.
7. **Fixed optical grid and model.** 400:700 nm at 1 nm; 3 illumination rows; ≤2 cameras; Hillman pathlength curve (`HemoCompute.m` L589-590). See section 8.

---

## 7. Current provenance

**Do outputs record which rig (or rig state) they were computed with? Essentially no.**

| Artifact | What it records about the rig | Evidence |
|---|---|---|
| `AcqInfos.mat` | `rigUUID`, `rigID` (binding only, no rig state/hash/version) | section 4 |
| `.dat` header (512 bytes, v1) | No rig field. Fields: magic, byte-order mark, version, flags, length, dtype, `frameRateHz`, `exposureMsec`, `channelName`, reserved bytes, axis sizes | `subFunc/datHeaderSchema.m` L51-58 |
| `HbO.dat`/`HbR.dat` (HemoCompute output) | Nothing about the rig. The header is built from the input channel header plus `outHeaderArgs` (L633-641). `opticalInfo.rigUUID/rigID` exist in memory (L1280-1281) but are not persisted | `HemoCompute.m` L589-641 |
| `dataHistory.mat` (PipelineManager) | Per output: `file`, `timestamp`, `pipelineVersion`, and per step `name`, `version`, `parameters`, `inputs(name,source)`. **No rig, no spectrum IDs, no rig `modifiedOn`.** HemoCompute's `FilterSet`/`Illuminations`/`HbT_concentration_uM`/`StO2perc` parameters *are* recorded | `classes/PipelineManager.m` L5770-5868 |
| `DataParams.mat → cameraCoregistration` | `resourceUUID`, `rigID`, a **copy of `tform`**, `sourceFile` (path), `sourceFileTimestamp`, `createdOn`, `appliedOn`, `appliedBy`, `method`, QC metrics. **No checksum and no rig UUID** | `subFunc/createDataParams.m` L201-233; set at `ImagesClassification.m` L246-279 |
| `ImagesClassification` second output `rigResolution` | Struct (`rigUUID`, `rigID`, `wasCreated`, `resolution`, coreg status) declared as a non-data pipeline output (L348-351). Whether it is persisted anywhere is UNKNOWN (not traced) | `ImagesClassification.m` L72-78, L348-351 |
| Rig side | `RigInfo.modifiedOn`; `.umit/logs/operations.log` (operation name + target UUID, no before/after); per-resource `checksum` for coregistration files | `UMITRigStore.m` L2706-2727, L2312 |
| Registration (not rig, but a precedent) | `DataParams.registration` stores the image reference's UUID and checksum used for the estimate | `createRegistrationTform.m` header, L8-12 |

**Can existing artifacts detect outputs computed with an outdated rig?** No automated mechanism exists. A human could manually compare an output's `dataHistory` timestamp against the rig's `modifiedOn` or `operations.log` (INFERRED; the log records operation name, not field-level changes).

**Pipeline step skipping.** `isHistoryMatch` (`PipelineManager.m` L13898-14017) reuses an existing output when step name, version, parameters (`isequal`), and input sources match. Rig state is not part of that match, so after an in-place rig edit (or rebinding) a re-run with `b_skipSteps` **would treat the old output as current**. (VERIFIED for the match rules in L13910-13917; the skip behaviour itself is documented in L4120-4128.)

---

## 8. Defaults and fallbacks

| # | Location | Default / fallback |
|---|---|---|
| D-1 | `UMITRigStore.getOrCreateDefaultRig` L1088-1100; `OiS200.json` | Empty/all-archived store → silently creates rig `OiS200` (2× D1024 cameras with `PF1024`, red `LED_632nm`, green `LED_521nm`, yellow `LED_593nm`) and makes it active |
| D-2 | `getOrCreateDefaultRig` L1070-1086 | Legacy `isDefault` marker promoted; or single active rig promoted to default |
| D-3 | `HemoCompute.m` L370-380, L386-409 | No binding → default/active rig (warning) |
| D-4 | `DatImageSource.m` L140-158 → `ensureDatasetRigAssociation` L480 | No binding → current Active rig, **written** into `AcqInfos.mat` on file open |
| D-5 | `UMITProjectStore.m` L1235, L1996 | `addSession`/migration → `ensureDatasetRigAssociation` (same Active fallback) |
| D-6 | `promptAddSessionDialog.m` L17-26 | Preview falls back to Active rig |
| D-7 | `DataViewer_Coreg2Cams › ensureRigContext` | No rig context → `getOrCreateDefaultRig()` |
| D-8 | `UMITRigStore.iUpgradeRigInfo` L2925-2927 | Missing `status` → `'active'` in memory |
| D-9 | `ImagesClassification.m` L488-493 | Missing camera model → `'D1024'` |
| D-10 | `importFromTif.m` L324-328; `ImagesClassification.m` L684-686 | Missing `CamIdx` → 1 |
| D-11 | `HemoCompute.m` L102-108 | `FilterSet='none'`, `Illuminations={'red','green','yellow'}`, `b_normalize=true`, `HbT_concentration_uM=100`, `StO2perc=60` |
| D-12 | `HemoCompute.m` L589-590, L599-613; `ioi_epsilon_pathlength.m` L59-64 | `'Hillman'` pathlength curve; 1 Hz smoothing and 1/120 Hz baseline in normalization; baseline `c_tot = HbT`; ±15 % concentration fit window |
| D-13 | `UMITRigStore.iResolveFilterSpectrum` L3210-3214 | Empty or `'none'` filter spectrum → flat `ones(1,301)` |
| D-14 | `getUMITRigSchema.m` L329, L336; `iValidateCameraRecords` L2999 | Fixed 400:700 nm grid, three canonical illuminations (`red`, `green`, `yellow`), ≤2 cameras |
| D-15 | `UMITRigStore.normalizeIlluminationName` L639-659 | `amber` → `yellow` |
| D-16 | `subFunc/resolveImportedChannelFallback.m` | Legacy folders: channel mapping inferred from filenames/headers |
| D-17 | `HemoCorrection.m` L83-86, L119-129 | Defaults `Red/Green/Amber = true`; channel files named `red`/`green`/`yellow` |
| D-18 | `GUI/resolveRigDisplayName.m`; `DataViewer.mlapp › resolveRigDisplayName` | Display name falls back to `rigID` when the store is unreadable |
| D-19 | `subFunc/createDataParams.m` L206-207; `loadDataParams.m` L134 | `cameraCoregistration.resourceUUID`/`rigID` default to `''` |
| D-20 | `UMITRigStore.iHydrateDefaultCameras/Illuminations` L3234-3278 | Built-in rig's display/manufacturer/model filled from spectrum metadata; `serialNumber` `''` |
| D-21 | `iAcquireWriteLock` L2585-2586 | `userName`/`hostName` default `'unknown'` in lock metadata |

---

## 9. Coupling with project management

Rig code itself has **no dependency on project code**: `RigManagement/*.m` only mentions `UMITProjectStore`/`getUMITProjectSchema` in comments (`UMITRigStore.m` L508; `getUMITRigSchema.m` L13, L79, L215). Every coupling below is **project/GUI code depending on rig code**, or data held in both.

| # | Coupling | Evidence | Technically necessary? |
|---|---|---|---|
| C-1 | Project session records store a copy of `rigUUID`/`rigID` (required session fields) | `getUMITProjectSchema.m` L252-253; comments L289-294 | **Incidental.** The authoritative copy is in `AcqInfos.mat`; the session value is a mirror. |
| C-2 | Project validation raises an **error** when a session has no rig (`session_rig_required`) or a half-filled pointer (`session_rig_incomplete`) | `UMITProjectStore.m` ~L5680-5700 | **Incidental** — a project policy, but it makes a "no-rig" session invalid at project level. |
| C-3 | `addSession` calls `ensureDatasetRigAssociation(saveFolder)`; caller-supplied `rigID`/`rigUUID` rejected | `UMITProjectStore.m` L1170-1176, L1235 | **Incidental** — side effect writes into `AcqInfos.mat` as part of project registration (stamps the Active rig). |
| C-4 | `updateSessionMetadata` rejects rig fields; `migrateSessionRigAssignment` and `iStampSessionRig` copy the dataset binding into the session record | `UMITProjectStore.m` L1914-1921, L1975-2000, L4795-4836 | **Incidental** (mirror maintenance) |
| C-5 | The mirror is only updated when the session record is blank. `assignDatasetRig` does not update it, and `migrateSessionRigAssignment` returns early when both fields are set (L1985-1988). DataViewer calls `refreshProjectContextForCurrentData` after reassigning (`SelectRigMenuSelected`), which reloads context but does not rewrite the project record | `UMITRigStore.m` L500-544; `UMITProjectStore.m` L1985-1988 | Consequence of C-1/C-4. A project-bound session's mirror can differ from the dataset binding after a reassignment (INFERRED: how `ProjectContext.rigUUID` is populated at L450 of the extracted DataViewer was not traced) |
| C-6 | DataViewer passes `ProjectContext.rigUUID/rigID` (the session mirror) to Coreg2Cams, but only in project-bound mode | `DataViewer.mlapp › OiSDUalCamCoregButtonPushed` | **Incidental** — the dataset's own binding would carry the same information |
| C-7 | `ProjectManagerTool` and `promptAddSessionDialog` display the dataset's rig and the Active rig as a preview | `promptAddSessionDialog.m` L15-27; `ProjectManagerTool.mlapp › (rig display)` | **Incidental** (display only) |
| C-8 | Save folder location: rig store is independent of project/save-folder location; `AcqInfos.mat` sits in the SaveFolder. No project path is used by rig code | `UMITRigStore.getRigsRoot` | Not coupled |
| C-9 | Rig checks during project setup: none by project code; project validation deliberately does not resolve the rig store ("to avoid coupling a project's own validation to another store's availability") | `UMITProjectStore.m` ~L5680-5686 | Deliberate decoupling |
| C-10 | Project's `UMITProjectBinding.umitlink` contains no rig data | section 4.6 | Not coupled |

---

## 10. Observations and risks

*(Suggestions are allowed here only.)*

**Latent bugs / inconsistencies**
- **O-1. Opening a dataset in DataViewer writes to it.** `DatImageSource` calls `ensureDatasetRigAssociation` (L151), which stamps the *current Active* rig into any legacy folder without a binding (D-4). For old data acquired on a different rig, this permanently records the wrong rig, silently, and is indistinguishable afterwards from a genuine binding. Suggestion: reads should never bootstrap a binding; use `readDatasetRigAssociation` and surface "unbound" in the UI.
- **O-2. `ImagesClassification` ignores the folder's existing rig** and always binds the machine default (L146-208). Re-importing into a folder bound to rig A binds it to the current default B. Only a warning when the default was auto-resolved.
- **O-3. `importFromTif` auto-creates `OiS200` without any warning** (L148), unlike `ImagesClassification` (L154-160). On a new machine every TIFF import silently invents a rig whose camera/LED spectra are those of one specific lab system.
- **O-4. Rig-store failures can block imports.** A malformed `defaultRig.mat` makes `listRigs` throw (L384), and a leftover `store.lock` has no clear API (L3353). Both stop `getOrCreateDefaultRig`, which every import uses.
- **O-5. "Default" and "Active" are conflated.** `setDefaultRig` = `activateRig`; `getDefaultRig` errors if the pointed rig is not active (L1010-1013); `HemoCompute.localResolveActiveRig` re-implements the fallback separately from `getOrCreateDefaultRig`. Two code paths decide "the current rig".
- **O-6. Rigs are mutable while referenced.** `setHardwareConfiguration`, `renameRigID`, `updateRigMetadata`, `updateFilterSet`, and in-place "Apply correction" change what an existing binding resolves to, with no version bump, no record of what changed, and no signal to outputs.
- **O-7. `renameRigID` breaks the `rigID` half of existing bindings.** `AcqInfos.mat` and `DataParams.cameraCoregistration.rigID` keep the old ID. UUID lookup still works, but `rigID`-only fallback (when UUID is blank) and `resolveRigDisplayName` fallbacks will fail. The stored `rigID` is never cross-checked against the rig resolved from the UUID: `ensureDatasetRigAssociation` opens by UUID (L469) and then compares the stored UUID with the opened rig's UUID (L484-487), i.e. a UUID against itself, so its "rigID and rigUUID identify different Rigs" error cannot fire when a UUID is present. A stale or wrong `rigID` therefore persists unnoticed.
- **O-8. Pipeline step skipping ignores rig state** (section 7). A HemoCompute output can be reused after the rig, its spectra, or its binding changed.
- **O-9. Reassignment leaves stale derived state** (section 4.4, C-5): coregistered `.dat` files, `DataParams.cameraCoregistration`, HbO/HbR, and the project session mirror all keep the old rig identity.
- **O-10. Coreg2Cams ignores the dataset binding** when not project-bound and uses/creates the machine default rig (section 5). A user can save a coregistration into a different rig than the one the data is bound to.
- **O-11. Spectra and filter sets are shared and referenced by ID only**, without a hash. User spectra are write-once by API, but `updateFilterSet` overwrites in place; built-in files ship with the toolbox and can change between versions. A binding therefore does not pin the optics.
- **O-12. Non-atomic `AcqInfos.mat` write.** `ImagesClassification.m` L692 uses plain `save([SaveFolder 'AcqInfos.mat'], …)` and then rewrites the file atomically with the rig fields (L208). Other writers use `saveMatAtomic`/`iSaveLoadedMatAtomic`.
- **O-13. `open` is O(number of rigs)**: it loads every `rig.mat` to match a UUID (L255-287), and several callers call it per access.
- **O-14. `GCaMP.json` uses the same spectrum (`FF01496LP`) for excitation and emission** (`resources/library/filterSets/GCaMP.json`); `none` and `jrGECO` have empty excitation. May be intentional (legacy mapping), but it is the filter set a typical GCaMP user will select. Worth confirming scientifically.
- **O-15. `iUpgradeRigInfo` defaults a missing `status` to `'active'`** (L2925-2927), which can make several legacy rigs "active" in memory; `validateLifecycle` then reports `multipleActiveRigs`; `create` refuses to run in that state (L216-220), `getActiveRig` throws, and `getOrCreateDefaultRig` throws (L1061-1063). `activateRig` is the only repair path (it demotes extra active rigs, L1481-1496).
- **O-16. Spectrum metadata ids/legacy aliases** exist in two places (`legacyMappings.json`, `migrateLegacyOpticalLibrary.m`); only the latter is exercised by tests as far as grepped. Low risk.

**Implications for the planned redesign (suggestions)**
1. *Run DataViewer without a rig.* Already mostly true: load paths catch rig errors (`DatImageSource` L150-157; `DataViewer` load alert). Remaining blockers are the side effects in O-1, C-2/C-3 (project requires a rig), and imports always resolving a rig (O-2, O-3; `ApplyRig=false` exists only for `ImagesClassification` and is not exposed to the pipeline).
2. *Snapshot binding.* The binding is the natural place for a snapshot: today it holds only identity. The data that actually affects computation is small: camera `index`/`spectrumID`, illumination `name`/`spectrumID`, the referenced spectra content (301 values each), the coregistration transform identity/checksum, and the filter set. A hash of that resolved optical configuration would be a stable "rig state" identifier.
3. *Provenance.* The `.dat` header has 423 reserved bytes (`datHeaderSchema.m` L58) but a header change is expensive (header spec). `dataHistory` step `parameters`/`inputs` and `DataParams.mat` are the more natural carriers; `isHistoryMatch` would need the rig-state ID as an input to treat stale outputs as non-reusable.
4. *Explicit no-rig binding.* Distinguish "never bound" (today: blank, and several paths fill it in silently) from "explicitly no rig". All the D-3…D-7 fallbacks would need to honor it.
5. *Immutability and export.* Content-addressed or versioned rig records would resolve O-6/O-7/O-11. The existing coregistration resource model (registry + checksum + active/archived folders) is a working precedent for versioned, checksummed items.
6. *Coupling.* C-1…C-7 could be removed without touching rig code. Project validation currently turns a missing rig into an error (C-2).

---

## 11. Open questions

1. **"Ledger" and "binding file".** Does the brief mean the rig store + `AcqInfos.mat` fields (what exists), or something planned? If a separate ledger is wanted, none exists today to preserve or migrate.
2. **Is `rigResolution`** (second output of `ImagesClassification`/`importFromTif`) persisted anywhere? Not traced (section 7).
3. **How `DataViewer.mlapp` populates `ProjectContext.rigUUID/rigID`** at the project-bound path (extracted doc L~450 builds `context` from project data). I did not trace whether it re-reads the dataset binding.
4. **Do tests cover** HemoCompute on an archived-rig-bound folder, a `renameRigID` after binding, or a corrupt `defaultRig.mat`? I grepped test names only; no match found.
5. **Intended semantics of `rigID` in bindings**: is it meant to be a human-readable alias only, or a fallback identity (current code uses it as a fallback when UUID is blank)?
6. **Is `GCaMP` (excitation = emission = `FF01496LP`) scientifically intended?** (O-14)
7. **Machine-default policy**: should imports on a store-less machine ever create `OiS200` automatically? (O-3, D-1)
8. **Concurrency**: the rig/store locks are folder-existence locks on a local path; multi-process or networked use of the rigs root was not examined.
9. **Unread areas**: `*_exported.m` files (off-limits by project rule), `EventsManager` (grep suggests no rig use; its many "rig" substrings are `trig…` identifiers), `docs/`, and the full `PipelineManager` apart from the provenance paths.
