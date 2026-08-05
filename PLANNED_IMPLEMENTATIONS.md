# Planned implementations

This file tracks fixes and improvements planned for the current development of the toolbox. It is not a roadmap or plan for a new toolbox version.

## High priority

### DataViewer pixel-selection event-patch refresh

- **Priority:** High
- **Problem:** Selecting a pixel calls `refreshEventPatches` through `refreshTemporalProfile`, even though an ordinary pixel change does not alter the event overlays. The current path deletes, recreates, and individually restacks every event patch on each click.
- **Profiling evidence:** With `fluo.dat` displayed at 512 x 512 pixels, three representative clicks spent 9.385 s in `refreshEventPatches` (98.7% of the 9.508 s click-callback time). `stackEventPatchesBottom` consumed 8.307 s across 252 `uistack` calls, while `DatImageSource.getPixelTrace` consumed 0.009 s total.
- **Planned implementation:** Keep event patches persistent during pixel selection. Remove the unconditional event-patch rebuild from the ordinary temporal-profile refresh path, and refresh or update patches only when event data, event selection, view mode, or relevant plot-axis limits change. Preserve the trace, SEM patch, time bar, legend, and graphics stacking behavior.
- **Validation:** Re-profile center and opposite-corner clicks on the same dataset, confirm that event patches remain visually correct, and verify that event-mode changes and event edits still rebuild overlays when required.
