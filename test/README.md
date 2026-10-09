# Test source and fixture policy

MATLAB test classes, runners, mocks, helpers, and fixtures under `test/` are tracked in version control.
Since the 2026-10-09 fixture trim the whole folder is about 280 MB, so the fixtures are tracked as
ordinary files (no Git LFS).

## Fixture directories

- `test/Analysis/TestingData_with_events/` (about 168 MB): headered `.dat` channels plus the full raw
  `img_00000.bin` and `ai_00000.bin`. It stays at full length on purpose.
  `testPixelValuesMatchHeaderlessReference` compares a SHA-256 of each channel built from the whole
  `img_00000.bin`, and `TestGetEvents` builds events from the whole `ai_00000.bin`.
- `test/Analysis/TestingData_retinotopy/` (about 64 MB): `fluo.dat`, `AcqInfos.mat`, and an `events.mat`
  whose `AnalogIN` is decimated 10x (`sr` is 1000 Hz). The raw `ai_00000.bin` is cut to one block.
- `test/Analysis/TestingData_speckle/` (about 46 MB): `speckle.dat`, `AcqInfos.mat`, a 3-frame
  `img_00000.bin`, and a one-block `ai_00000.bin`.

The speckle and retinotopy raw files are structurally valid but short: no test reads them in full, and
tests copy only the `.dat`, `AcqInfos.mat`, and `events.mat`. `test/Tools/trimLargeFixtures.m` documents
exactly what was trimmed and why. The full recordings were removed on 2026-10-09 and are not kept in the
repository.

Keep fixtures small. Prefer synthetic data built in a temporary folder
(see `TestSpeckleMappingNaNDenominator`) over adding recordings. No single file may exceed 100 MB, the
GitHub hard limit; the largest is `TestingData_with_events/img_00000.bin` at about 76 MB.

## Local diagnostics

`failed_tests.txt` files are local diagnostics, not current validation evidence, and are Git-ignored.
Use `_tools/validateTaskChanges.m` and its current report for test status.
