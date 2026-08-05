# Documentation rebuild feedback
- In case of conflict, comments in the sections supersede general comments.
## General comments:
- Avoid coding jargon like "schema-validated file".
- Remove the following sections from all Getting Started pages:
    - troubleshooting
    - Expected result
    - what is saved: replace this one by a note box.
- Rewrite "Prerequisites" header. Maybe just remove the header and replace by a sentence like:
" To execute this function, first be sure to load a .dat file...".
- Same fof "Purpose" header. Remove the header and rewrite the sentence. Example: replace "Build and run a small, reproducible sequence of analysis functions on the active recording." with "The Pipeline Manager tool is responsible to Build and run a small, reproducible sequence of analysis functions on the active recording."
- When describing the operation sequences, add inline icons when referring to existing tools that are shown in DataViewer's tabs. Example: Open the Pipeline tab <INSERT ICON HERE> and choose Pipeline Manager.
- For inline code, ensure font is Monospace or other that is more distinct from the text font.
- Remove the header ### Purpose from all pages in "Getting Started". Keep only the section's body text.
- Put all apps and functions names in the text in Italic.
- Update systems requirements to allow Matlab's version R2024b and onwards. Other parameters are the same.
- A brief explanation of what a .dat file and a .umt file is needed, mainly in DataViewer pages. An introduction to this two-file framework is needed. The main idea behind having .dat and .umt stems from the fact that now, all .dat files in the SaveFolder have the same dimensions and these are managed by the AcqInfos.mat file. Therefore, .umt files allows to store other kinds of data that diverge from AcqInfos.mat and contain the whole meta data associated with it. In a typical imaging data processing workflow, .umt files will normally contain processed/reduced imaging data (with way less volume than the raw .dat files) and represent intermediate or final steps in a data processing pipeline. Some of those can be opened by DataViewer (image kind).
- Some tables frequently emphasizes if an operation changes the data or is persistent or temporary. This is not needed everywhere. At the Documentation of DataViewer, add a general explanation about what is kept and what is saved. 
- To expand each section, read the corresponding older documentation for inspiration. Do not take for granted that the old documentation is still relevant. Cross-validate the statements in the old documentation with the current implementation. In doubt, do not use the old documentation to create text in the new one. 
## Home
ok
## Getting Started
ok
## Introduction to DataViewer
### Purpose
- Replace:"DataViewer is the main workspace for opening mesoscale imaging data, moving through frames and events, inspecting pixel and ROI signals, and sending the current recording to processing and project tools." With: "DataViewer is an app dedicated to exploring mesoscale imaging data. The app allows one to perform basic data exploration by moving through frames and events. One can also inspect the temporal profile of individual pixels or ROIs."
-Remove ### Prerequisites section
-Remove ### Expected result section
-Put ### What is saved section under a note bubble.
## Opening data
-Remove ### Prerequisites section
### Import a recording
- Add a comment that bin files correspond to the data created using the OiS200 LightTrack widefield imaging system (link to https://labeotech.com/product/modular-optical-imaging-system/)
- Tif files are generic files for other rigs
- Remove ### Expected result section
- Put ### What is saved section under a note bubble.
## Exploring imaging data
- Replace ### Prerequisites header by "Basic operations"

## Working with ROIs
- The ROI creation and edition workflow lacks details. For example, to confirm ROI creation double-click inside it. 
- Remove ### Prerequisites section.
- Add Section header to the 5 operations.
- The section about Allen Atlas procedure is too detailed to be here. Just point to the ROI management documentation for details. Make a high-level description of what the tool does.
- The description of ROI creation by thresholding is lacking here. Same comment here as for the Allen Atlas detail level.

## Exploring temporal traces
- This section lacks images to illustrate ROI table creation and ROI traces

## Basic processing workflows
- Describe briefly the layout of the PipelineManager tool (Available functions on the left, 
Graphical representation of the processing pipeline in the center, Step Params on the right)

## Documentation
ok
## DataViewer
ok
## DataViewer overview
- replace "opens continuous .dat recordings and image-kind .umt data" by "opens image time series"
- replace ### Supported sources by ### Supported files. Add a link for more information about each data file in the table
- Remove dialog box "Optional metadata degrades safely".
- Rewrite text of "Persistence at a glance". Jargon is too technical. Convey this idea: "DataViewer is conceived to perform basic data exploration. So the user may not need to save everythin. It automatically saves only necessary metadata for display, pixel coordinates and mask. Processed data, newly created ROIs etc are not automatically saved.

## User interface
- Page lacks illustrations. Add at least one screenshot.
- Remove "Utility dialogs" section. Instead, expand the "Tabs" section to include a brief description of each element of the tabs
- Add links to the corresponding documentation in the Tabs table (controls) when the dedicated documentation pages exist.
- Remove "Main work areas" section. Instead convey the same info from the table as text under a screenshot.

## Data loading and navigation
- Loading sequence section: This should be more like a step by step tutorial for the user rather than a description of what is done unde the hood. Create a sequence like: "To load an image file in dataviewer, go to File > open and select either a .dat file or a .umt file..."
- Expand "Navigation controls" explaining that Event mode is available only when events were detected at data import and an events.mat file exists in the Save Folder (refer to EventsManager documentation for details).
- Change scope of this page to data loading and exploration. Try not to be redundant with corresponding section in "Getting started". See related comment in "Image visualization" below.

## Image visualization
- Merge this page with "Data loading and navigation". Change the name of the latter to "Data loading and basic exploration"
- Change "Spatial units" to "Image Calibration". Provide a more detailed explanation of the calibration options. Add a screenshot of the image calibration dialog box.
- "Cache behavior": Explain that the cached image region will be shown inside a dotted white rectangle. By clicking outside the rectangle the cached area will be updated automatically. 
## ROI management
- "Create ROIs" section lacks screenshots. 
- A section describing the procedure to draw rectangles, ellipses and polygons is lacking.
- In "Allen atlas selection and placement", instead of  "Use a module context menu", say "right-click on a module in the Tree".
- In the "Scientific limitation" box, add the following caveat note: "The position of the Bregma and Lambda is provided here as a rough estimate. The data from the Allen Mouse Brain Atlas does not provide any anatomical landmark coordinates from the mouse skull. For more info on this, see this discussion (link to https://community.brain-map.org/t/why-doesnt-the-3d-mouse-brain-atlas-have-bregma-coordinates/158) from the Allen Brain Map Community Forum. 
- "ROI table and editing" This section needs more detailed explanations about: color changing, name changing, multi-selection and ROI boolean operations.
## Temporal visualization
- Merge this page into "Data loading and navigation"
- Condition/Repetition deletion is not for display-only. It updates a list of ignored events in the events.mat file that will continue to be excluded even during data processing. Be sure to convey this message.
## Processing and export
- In the first paragraph, explain that the data processing is managed by the Pipeline Manager tool. Link to the related tool.
- The concept of virtual file-source node needs a brief explanation or change the term to a clearer one.
- Create a dedicated section to "Data History" explain what the Data History shows as information. 
- "Temporary output lifecycle". Explain that, by default, image results are stored as temporary. If the user doesnt choose to save, it will be deleted from the SaveFolder.
- "Other output commands": Change the section header to something like "Data Export". 
- Remove "Alignment is different" box. A dedicated page is needed for Image Reference creation/management and Image Alignment. It should be placed after "Project and session integration" page.
## Project and session integration
- This section is too technical and explains too much about the details of the project management.
- This section should explain the high-level reasons why we would need a "Project". THe main reason here is to track multiple recordings from a single mouse or from a group of mice throughout the evolution of a research project. The project management tool allows the user to have a structured organization of his project without the need to have specific Save Folder file naming conventions to organize their imaging data. A project also allows to  structure the creation of Image References for a mouse. This is particularly relevant to longitudinal studies where the same animal is recorded several times. Frequently, one would want to align all recordings of an animal so all the data falls inside the same spatial coordinates. Image reference creation and management and Image Alignment tools do just that. Project creation is optional but highly recommended if the experiment starts to scale up from single mouse pilot studies towards mouse cohorts and longitudinal studies.
- Just explain briefly that the Project management doesnt rely on Folder Names or forces the user to save the data in a particular folder or forces a folder tree structure of any kind. It uses the .umitlink in the SaveFolder to locate the related metadata in the project.
- Explain that the Project info is "centralized", meaning that it is saved to a root folder (show command to see it) in the user's root. The project is structured as follows: a Project contains one or more Subjects. Each subject contains one or more recording sessions. Each of those will have specific meta data associated with it. 
- A project can be set and managed using the dedicated tools (add links);
- Remove "Ownership boundary" section.

## Related Apps and Tools
ok
## Pipeline Manager
- This page needs to be carefully expanded to provide more details on the GUI operation as well as on the basic and advanced option settings.
- First paragraph: Explain that the Pipeline Manager allows one to build a DAG-like sequence of analysis functions. The existing functions cover most cases for the analysis of widefield imaging datasets. For custom functions, refer to the related advanced User resource. Explain that the manager exist to ensure the tracking of the operations of the data as well as the optimization of RAM usage. It also provides a graphical representation of the pipeline which facilitates the visualization of the flow of data between operations (steps).
- Add a syntax section to show how to open outside DataViewer as a standalone app.
- "Interface": add a brief text description of the app layout: the analys functions are listed on the left, a graphical representation of the pipeline is listed on the middle, the step summary on the right and main action buttons on the top. 
- Rewrite: "Selects leaf persistence. DataViewer locks exploratory runs to viewer-managed temporary output where applicable."  This is not clear. Put it in more lay terms like: Control automatic saving of last steps. Disabled for DataViewer.
- "Basic workflow": This is not very user-friendly. Instead use a tutorial-like approach: "1- select a function from the list and click add... 2-Depending on the function input needs, a dialog box will appear prompting to selec the input to the function... 3- in the graph, right click on the step (function) to access "Edit Parameters" ...". This section lacks illustrations.


- "Outputs and DataViewer refresh": The text is not clear. This explanation seems to refer to DataViewer usage only. In this case, explain that when called from DataViewer, image-compatible results are displayed, otherwise the data is saved as .umt files in the SaveFolder when.
- Explain briefly the content and structure of files that are controlled by the Pipeline manager (dataHistory.mat and pipeLog.mat). Why they exist and what info they contain.
- Merge "Files and state" to "Outputs and DataViewer refresh", then change the latter to reflect the content (output files and control files).
- Specific sections are needed with more detailed workflow and purpose explanation:
    - Graph interactions: step selection, parameters editing, step reconnection and deletion etc.
    - Pipeline Script Generation
    - Data Folder Selection
    - Execution Settings
    - Advanced RAM options
    - Tools menu: explain individually each submenu option.
    - Adding steps as new branch: explain that one can expand the pipeline to process multiple files with a single pipeline build by creating parallel branches.
- "Function warnings still apply." box: The text is not clear. Give examples like " Data import functions overwrites the folder content so be careful when executing these functions". Advise to read the documentation. 
## Events Manager
- First paragraph: state that the app is dedicated to the visualization and detection of events stored in the analog inputs from LabeoTech Imaging systems (ONLY, for now). It enables users to easily set the parameters for event detection either from internally generated events from the imaging platforms or from external sources such as TTL signals or photodiodes.
- Use the image "events-manager_too_edgeType.png" from the corresponding old documentation to illustrate the trigger type.
- Detection parameters need further explanation (see old documentation).
- Explain how trials are splitted. The current approach is that one sets the Baseline period (i.e. a period before the event onset) while the total trial length is the sum of the baseline period and the time between the event onset and the frame before the next event onset. 
- A more detailed explanation about the condition file selection is needed.
- Underscore the idea that the goal of the app is to ultimately create an events.mat file in the SaveFolder. This file will be used by DataViewer to set the events view of .dat files as well as the analysis functions will use it as a source of information to split the data by events.
## Image Reference Manager
- Replace current screenshot with one that contains image references. I will provide info about which file to use.
- "Basic workflow": add that the image can be edited (affine transformation, confirm) as an optional step.
- "Identity and history" box: add "only" after "manager". Put the "Do not" in bold.

## Image Alignment Tool
- Replace the screenshot with a better one with a full set of Image reference + moving
- "Automatic registration": link to a matlab documentation of the main function used to do the registration
- "Reset Fine Tuning /Reset Registration: "Removes manual adjustments or the full candidate." is not clear. Rewrite to something like:" Resets changes".
- "Apply Registration to Folder": Explanation is too technical. Say that the geometric transformation will be applied to the data in the Save Folder. point to the Note box below.
- Put the two comments below in a dedicated section:
    - Add an explanation about which files in the folder are affected by the transform and which are not.
    - Add an explanation about where the transform (i.e. tform object) is stored and how the user can know if the folder was already geometrically transformed or not.
- Add a section explaining the QC (quality control) metrics applied in the coregistration and where those are stored.

## Dual-Camera Coregistration
- Replace first paragraph with "This tool is dedicated to perform the coregistration of the images from Dual-Camera Imaging systems from LabeoTech. The goal of this tool is to ensure that the output of both cameras are aligned."
- Add a dialog box stating that this coregistration is normally a one-time calibration procedure for dual-camera recording systems. It is advisable to perform this operation when the hardware is moved or its optics are changed.
- "Important source rule" : remove "The Save Folder passed by DataViewer is context only." Underscore the fact that an unregistered set of imported .dat files is needed for this.
- I noticed that there is no way to perform a second calibration procedure without manually deleting the calibration file in the root folder. This is not acceptable. Apply the above changes and flag this doc update as deferred and depending on GUI refactoring.
- Add an explanation about which files in the folder are affected by the transform and which are not.
- Add an explanation about where the transform (i.e. tform object) is stored and how the user can know if the folder was already geometrically transformed or not.
- Add a section explaining the QC (quality control) metrics applied in the coregistration and where those are stored.
## Project management and folder binding
- Replace jargon " project store" with "project info" or similar
- "Project Manager" section: I cant find the optiosn "Set Active/ Archive / Restore Resource". Check if these options really exist, if so, improve explanation. If not, remove from table.
- "Unbinding, recovery and repair": Rewrite in a more accessible vocabulary.
- "DataViewer refresh": remove "Project metadata never absorbs the bulk imaging files."
- "UUIDs are identity; paths are locators.": change box title. It is too technical and descriptive. Convey call for action instead.


## Functions
ok
## Image registration
ok
## applyRegistrationTformOnFolder
- State which file types are modified by the transform and which are not by thre transform.
- NOTE TO REFACTO THE FUNCTION ITSELF: Force the DataParamsFile to be DataParams.mat only. Remove this from Parameter. Set this operation as planned and update this documentation when it is implemented.
## createRegistrationTform
- Create detailed explanation of quality control metrics and how to interpret them.
## genAmplitudeMaps
- No need for the Low-RAM .dat behavior box.
## genImageTimeSeriesUMT
- "Validation includes event rules." box: The text is not clear. "Keep labels and any later event metadata aligned with the entry value.": what does it means exactly?
## Advanced User Resources
ok
## Events
- Add an asterisk to "AnalogIN" field. State that it is specific for the Labeotech imaging systems.
- State that this file is created and managed by the class EventsManager. The app EventsManager tool is a gui wrapper of the class and provides easy control of the options for events management (link to documentation).
- "Current events.mat variables" table: add the following information about each field: data class, dimensions.
- Document that external signals (outside LabeoTech systems) can be detected using the EventsManager class direclty in matlab. Currently this option is not supported in the corresponding app. Show an example of code using detection of external signal.
- "How DataViewer uses events": State that Events Manager stores flags the 'deleted' conditions and repetitions. This is not a destructive operation. Instead, the trials set as "deleted" are simply ignored during visualization in DataViewer and during any analysis function that uses the events.mat file.
- Provide more details (not much) on how event information is stored in umt files or link to advanced user resource page. Create a small table with the fields and their meaning, data class, dimensions.
- "Recovery" section: Refrain from explaining what happens with the Events Manager tool. This page focus on the Events file. Put this in a dialog box: "Replace an existing events.mat only through the tool or validated code, and keep backups for irreplaceable event annotations.".
## Save Folder metadata files
- "A Save Folder keeps recording data and recording-local state together. Centralized project identity and managed resources remain outside it.": Add that although the project resources live outside, a link file exists. 
- Underscore that all .dat files are SaveFolder-bound where their dimensions are the same and stored in AcqInfos.mat file. 
- Underscore that each meta data file needs to be set through their own dedicated class/app. Do not attempt to hand-edit if the user doesnt know what he is doing.
- Add important note box: Removing AcqInfos.mat will render .dat files un-readable.
- Consider the first table in the page as a summary. Add a detailed explanation of the structure of each meta data file below.
- Remove *.dat and *.umt files from the table. They have their own description pages.
- Remove *.roi from this page and create a dedicated description page under Advanced User Resources.
- Remove PMTMP_* from the table.
- "Centralized versus local" section: rewrite the text and avoid terms as "authoritative" and "store". Make a more hihg-level explanation so the user knows the overall framework of the toolbox without too mych details about UUIDs and others.
## Organization of .dat files
- Add this information in the first paragraph : "Imaging datasets are stored in binary files with the extension .dat. Each .dat file contains only a single channel. For example, in a recording where a fluorescence channel was acquired with red and green illuminations, each channel is stored separately in a dedicated .dat file. These files are created using one of the available Data Import functions. The binary files store only imaging data as 4-byte (32-bit) floating-point values (single precision)."
- "Metadata and file naming" section": Remove "A multi-camera Save Folder may contain camera-specific source names and several derived channel files, each reconstructed with matching metadata."
## Organization of .umt files
- Explain why .umt files exist. 
- "Entries": State that each .umt file may contain different entries with distinct dimensions. See 'entries' as 'measures' in the context of a processed dataset that derives different measurements from it. 
- "Create, append and validate": State that normally this is done by PipelineManager and the analysis functions. How the example codes for advanced use.
- "Identity boundaries": put this into a note box.
- "Backward compatibility and editing": state that data generated using older versions of the toolbox, mainly processed data saved as .mat files, may not be compatible with the .umt structure.

## Creating analysis functions for PipelineManager
ok
## Archived (previous releases)
ok
## Previous Getting Started
ok
## Previous Apps
ok
## Previous Tools
ok
## Previous Advanced User Resources
ok
