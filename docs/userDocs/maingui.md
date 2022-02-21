### The main user interface
___

The *Umitoolbox* main graphical user interface (GUI) allows one to create and manage a project and its datasets, configure and run analysis pipelines and launch visualization apps.

##### Syntax
___

`Umitoolbox`: Opens an empty interface. Create or load a [project file](/creating_a_new_project.md).

`Umitoolbox(ProjFile)`: Opens the interface and automatically loads the project stored in `ProjFile` .mat file. The .mat file `ProjFile` is a string with the full path to the project file.

##### Sections
___   

* [The interface](#the-interface)
* [Advanced Search](#advanced-Search)

#### The interface
___

<p align="center">
  <img alt="guiMainComponents" src="../../assets/img/mainGUI_mainComp.png"/> <br>
  <em>Main graphical user interface.</em>
</p><br>

##### 1. Menu bar
**Project**
* ***New:*** Creates a new *project file*. Here, the user will be asked to create a project name and to select the save and the raw data folders as well as to select the protocol function associated with the raw folder structure. For more info, see [*Creating a new project*](/creating_a_new_project.md).
* ***Open:*** Opens dialog to load previously saved *project file*.
* ***Save:*** Saves the currently opened *project file*.
* ***Save as...:*** Opens dialog to save currently opened *project file*.   

**Utilities**
* **Project >>**
    * ***Look for new data:*** Scans the raw directory for new/deleted recording folders.
    * ***Update Project:*** Updates the save folder subdirectories when new recordings are added to the raw folder.
    * ***View LogBook:*** Shows a table (data from *LogBook.mat* file in the save directory) with all operations performed in the data so far.
    * ***Reset LogBook:*** Erases the LogBook table.
    * ***Restore garbage list:*** Erases all objects listed in the *garbage list*. For more info, see [Managing experiment info](/ht_manage_exp.md).

##### 2. Search tool panel
Provides options for filtering the data shown in the *object tree*.

* ***Advanced Search:*** Launches an interface to filter the project items (Subjects/Acquisitions/Modalities). For more info see [Advanced Search](#advanced-search).
* ***Reset Search:*** Resets the *object tree* to show all items of the project.

##### 3. Tabs' options panel
This panel displays options for each tab. The content of the options panel changes depending on which tab is selected. For details on each tab's options, read the respective tab documentation in [section 5](#tabs).

##### 4. Project's object tree
Shows a tree containing the file structure of the project's save folder. The tree behaviour changes depending on the selected tab. *Dashboard* and *Visualization* tabs allow single-item selection while the *Pipeline control panel* tab allows multi-selection. The object tree is disables for the *Analysis* tab.

##### 5. Tabs
The content of each tab.
Here are more detailed information on the components of each tab:
* [Dashboard](/maingui_dashboard.md) (Experiment Management)
* [Pipeline control panel_LINKNOTFUNCTIONAL](/maingui_pipeline.md) (Data processing)
* [Visualization_LINKNOTFUNCTIONAL](/maingui_visualization.md) (Data visualization)
* [Analysis_LINKNOTFUNCTIONAL](/maingui_analysis.md) (Results visualisation and analysis)

#### Advanced search
___

Content available soon!

\
[**<< Home**](../../index.md)                                                                   
