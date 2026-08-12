# Static analysis report: `Analysis/`

Date: 2026-08-11  
Scope: all 28 `.m` files under `Analysis/`  
Method: MATLAB Code Analyzer plus read-only source inspection. No unit tests, runtime execution, or code edits were performed.

## Executive summary

- No MATLAB Code Analyzer error-severity findings were reported.
- No definite syntax failure was identified by the static analyzer.
- Three conditional correctness risks merit focused review before relying on unusual event metadata or boundary-shaped recordings.
- The remaining analyzer findings are warnings about unused arguments, obsolete suppression comments, deprecated date APIs, or performance—not confirmed functional defects.

## Potential correctness issues

### 1. Event IDs are used as direct indices into event-name lists

**Severity: medium; conditional on event IDs not being contiguous 1-based indices**

- `Analysis/Imaging/genAmplitudeMaps.m:190`
- `Analysis/Imaging/genAmplitudeMaps.m:350`

Both paths construct names with `eventNameList(eventIDs)`. The surrounding validation requires positive integer IDs but does not establish that IDs are contiguous or bounded by the name-list length. If an `events.mat` file contains IDs such as `[10, 20]`, this can either throw an indexing error or associate names incorrectly. The computation itself groups by the IDs, so name lookup should be based on an explicit ID-to-name mapping rather than assuming IDs are list positions.

### 2. Retinotopy averaging mode does not validate matching or sufficient onset/offset events

**Severity: medium-high; input-dependent runtime failure**

- `Analysis/Retinotopy/genRetinotopyMaps.m:220-226`
- `Analysis/Retinotopy/genRetinotopyMaps.m:240-244`
- Low-RAM equivalent: `Analysis/Retinotopy/genRetinotopyMaps.m:279-287`, `:323-325`

The code indexes `indxOn(2:end)`, `indxOff(1:end-1)`, and `indxOff(ii)` without first asserting that onset and offset lists are nonempty, equal in length, and contain enough transitions. With a missing offset, an unmatched transition, or only one onset, the function can produce empty/`NaN` lengths or index beyond the available offset list. In average-movie mode, the low-RAM path additionally clamps the start frame to 1 but later indexes `trialData(:,:,1:bsln_len)`; recordings whose first baseline extends before frame 1 can therefore fail with an out-of-bounds index.

### 3. Response-feature extraction permits a baseline that consumes the whole recording

**Severity: medium; boundary input can produce an empty response window**

- `Analysis/Stats_Functions/calculateResponseFeatures.m:114-125`

`baselineFrames` is clamped to `nT`, but no check requires at least one frame after the baseline. If `baselinePeriod * FrameRateHz >= nT`, the default response window becomes `baselineFrames + 1:nT`, which is empty. Subsequent `max`, latency, and feature calculations then operate on empty arrays and can return invalid outputs or fail. The function should reject this input explicitly or define a documented empty-response policy.

## MATLAB Code Analyzer findings

These are reported as potential maintenance/performance issues, not confirmed runtime defects.

### Warnings

- `Analysis/Stats_Functions/calculateResponseFeatures.m:145` — use `strcmpi` instead of applying `strcmp` to normalized text. Functional behavior is currently equivalent because the value was lowercased earlier.
- `Analysis/Retinotopy/genRetinotopyMaps.m:205` and `:257` — `SaveFolder` is unused in the two execution helpers. Confirm whether it is intentionally retained for the shared function signature or should be replaced with `~`.

### Informational performance/modernization findings

- `Analysis/Imaging/split_data_by_event.m:226` — scalar detection can use `isscalar`.
- `Analysis/Imaging/run_HemoCorrection.m:105` — scalar-length comparison can use `isscalar`.
- `Analysis/Imaging/genAmplitudeMaps.m:303` and `:467` — `find` results are used only for indexing; logical indexing may be faster.
- `Analysis/Filters/apply_aggregate_function.m:320` — same `find`/logical-indexing advisory.
- `Analysis/DataImport/getEvents.m:360,370` — cell-content extraction can use direct indexing instead of nested cell indexing.
- `Analysis/Filters/normalizeZScore.m:139` — nested cell indexing can be simplified for performance.
- `Analysis/ImageRegistration/createRegistrationTform.m:121` — `ismatrix` is preferred when checking for a matrix.
- `Analysis/ImageRegistration/createRegistrationTform.m:252,304,323` and `Analysis/ImageRegistration/applyRegistrationTformOnFolder.m:191` — `now`/`datestr` are deprecated in favor of `datetime`-based APIs.

### Stale suppression comments

The analyzer reports obsolete `#ok`/suppression directives in these files. They do not indicate a current code defect, but they can hide future diagnostics and should be cleaned up during maintenance:

- `Analysis/Stats_Functions/getDataFromROI.m:59`
- `Analysis/Stats_Functions/genImageTimeSeriesUMT.m:22`
- `Analysis/Stats_Functions/genCorrelationMatrix.m:153`
- `Analysis/Stats_Functions/calculateResponseFeatures.m:38`
- `Analysis/Retinotopy/genVSM.m:45`
- `Analysis/Retinotopy/genRetinotopyMaps.m:39,275`
- `Analysis/Imaging/split_data_by_event.m:32`
- `Analysis/Imaging/run_SpeckleMapping.m:107,114`
- `Analysis/Imaging/run_HemoCorrection.m:94`
- `Analysis/Imaging/run_Ana_Speckle.m:31`
- `Analysis/Imaging/genAmplitudeMaps.m:37,81,294`
- `Analysis/Imaging/apply_detrend.m:51,369,374`
- `Analysis/ImageRegistration/createRegistrationTform.m:516,520,524,528`
- `Analysis/ImageRegistration/applyRegistrationTformOnFolder.m:151,158`
- `Analysis/Filters/spatialGaussFilt.m:314,319`
- `Analysis/Filters/normalizeZScore.m:52,230,235`
- `Analysis/Filters/normalizeBSLN.m:477`
- `Analysis/Filters/GSR.m:334,342`
- `Analysis/Filters/apply_aggregate_function.m:455`
- `Analysis/DataExport/run_ConvertToTiff.m:26,201`

## Files with no Code Analyzer findings

`run_HemoCompute.m`, `funcTemplateUMT.m`, `funcTemplateFileManifest.m`, `funcTemplate.m`, `normalizeLPF.m`, `run_ImagesClassification.m`, and `importFromTif.m` returned no Code Analyzer issues in this pass. This is not equivalent to proof of correctness; it only means the built-in static analyzer reported no diagnostics for those files.

## Limitations

This report does not execute functions, load fixtures, inspect runtime types, resolve every toolbox dependency, or verify numerical results. The three correctness items above should be confirmed with focused tests or representative data in a separate task.
