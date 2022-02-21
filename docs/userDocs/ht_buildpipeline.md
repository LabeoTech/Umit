### How to build an analysis pipeline

One can build an analysis pipeline through the main GUI [Pipeline control panel](/pipeline_tab.md) tab or using the standalone version of the [DataViewer](/dataviewer.md) app.
##### The Pipeline Bulder interface
___   

![pipelineBuilderGUI](/assets/img/pipelinebuilder_mainFig.png)

##### 1. Function List
List of the available analysis functions separated by category. These functions are located in ```./umIToolbox/Analysis```.
Functions that have optional parameters are indicated by an exclamation point ```(!)```.

##### 2. Add function options
* ***Add function to:*** Dropdown menu to choose which object from the *Protocol* will be used in the selected function. Not applicable to standalone version of ```DataViewer```.
* **Add button:** Adds the selected function from the *Available functions* list to the *Selected functions* list.

##### 3. Pipeline
List of functions (steps) from the pipeline. The pipeline will run as the order shown in the list (from top to bottom).

##### 4. Step controls
* ***Up/down arrows:*** Move the selected step up/down in the list
* ***Set Options:*** Opens a dialog box containing extra parameters for the selected function. Check the documentation of each function in the [analysis functions](/index.md/#analysisfunctions) section for details on optional parameters.
* ***Remove item:*** Removes the selected function from the pipeline.
* ***Save step:*** Saves the selected step to a .dat file. Opens a dialog box to type the name of the .DAT file. Some functions do not provide imaging data as outputs and will not allow saving the step. Please, refer to the [analysis functions](/index.md/#analysisfunctions) section for details.

##### Running an analysis pipeline from *DataViewer* standalone app
___   

To execute an analysis pipeline in the standalone version of the ```DataViewer``` app:
1. Open the imaging file to be processed.
2. Click on *Utilities >> Data Processing* to open the pipeline builder interface.
3. In the pipeline builder interface, select the analysis functions to run on the file.
4. Click on *Create Pipeline*. 
5. A window will appear containing a summary of the pipeline steps. Click on *Run* to run the pipeline.
6. Once finished, the final processed file will be displayed. If the last step of the pipeline was not saved during the pipeline creation in the Pipeline Builder interface, the file will be saved with a temporary file name with prefix ```tmpFile_```. You can save the processed imaging by clicking on *File >> Save as...*. 


> Note: All temporary files (```tmpFile_xxxx.dat```) contained in the save directory will be deleted by closing the app. Be sure to save the displayed imaging data before closing.


##### Running an analysis pipeline from umIT's main interface
___  

To execute an analysis pipeline from the main GUI of **umIT**:

1. In the *pipeline control panel* tab, select the objects to be processed in the folder tree.
2. Save selection by clicking on *Select* button.
3. Click on *Configure New Pipeline* to open the pipeline builder interface.
2. In the pipeline builder interface, select the analysis functions to run on the file.
3. Click on *Create Pipeline*. 
4. In the *pipeline control panel* tab, a summary of the pipeline will be displayed on the *Pipeline status* panel.
5. Click on *Run Pipeline* to execute the pipeline over the selected objects.
6. Once finished, a table containing a summary of the pipeline will be displayed in the *Pipeline Summary* panel.














[**<< Home**](/index.md)
