# Planned implementations

This file tracks fixes and improvements planned for the current development of the toolbox. It is not a roadmap or plan for a new toolbox version.

## High priority

## Medium priority

### Standardize the registration metadata filename

- **Priority:** Medium
- **Problem:** `applyRegistrationTformOnFolder` exposes a `DataParamsFile` option even though Save Folder ownership requires registration state to be stored in the standard `DataParams.mat` file.
- **Planned implementation:** Make `DataParams.mat` the fixed registration-metadata target. Remove the `DataParamsFile` name-value argument from the input parser, help text, and PipelineManager metadata, then update the function documentation to match.
- **Validation:** Update and run the focused registration-function tests, verify PipelineManager discovery and option rendering, and confirm that folder registration still updates `DataParams.mat` while rejecting accidental reapplication.
