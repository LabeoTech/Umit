### How to build an analysis pipeline

One can build an analysis pipeline through the main GUI [Pipeline control panel](/pipeline_tab.md) tab or using the standalone version of the [DataViewer](/dataviewer.md) app.
##### The interface
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

##### Building a pipeline in *DataViewer* standalone app
___   









[**<< Home**](/index.md)
