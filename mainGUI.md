### The main user interface

The *umIToolbox Main* graphical user interface allows one to create and manage a project and its datasets, configure and run analysis pipelines and launch visualization apps. 
#### Main elements
The interface is constituted of a menu bar and four tabs. Here are the main components of the interface:
![fig1](/assets/img/umIT_mainGUI_mainComp.png)
###### Fig 1. Main GUI components.   

##### 1. Menu Bar
**Project**
* ***New:*** Creates a new *project file*. Here, the user will be asked to create a project name and to select the save and the raw data folders as well as to select the protocol function associated with the raw folder structure. For more info, see [How to create a new project](/ht_create_new_project.md).
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

##### 2. Tabs
Each tab deals with a one specific role of the toolbox (i.e., experiment management, data analysis, visualisation and statistical analysis). See links in section 5 ([Tab content](#tab-content)) for details.
##### 3. Search Panel
Provides options for filtering the data shown in the *folder tree*.
!Currently, each tab as its own *search panel* and *folder tree*.
* ***Advanced Search:*** Launches an interface to filter the project items (Subjects/Acquistions/Modalities). For more info see [Advanced Search](/ht_use_advanced_search.md).
* ***Reset Search:*** Resets the *folder tree* to show all project items.

##### 4. Folder Tree
Shows a tree contaning the file structure of the project's save folder.
##### 5. Tab content
The content of each tab. 
Here are more detailed information on the components of each tab:
* [Dashboard](/dashboard_tab.md) (Experiment Management)
* [Pipeline control panel](/pipeline_tab.md) (Data processing)
* [Visualization](/visualization_tab.md) (Data visualization)
* [Analysis](/analysis_tab.md) (Results visualisation and analysis)

[**<< Home**](/index.md)                                                                       [**Next item >>**](/file.md)