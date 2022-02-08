# Welcome to the umIT (*u*niversal *m*esoscale *I*maging *T*oolbox) Wiki page!  
Recent advances in genetic tools, such as fast-responding calcium indicators (e.g., GCaMP6f) as well as the increasingly accessible high-speed, high-resolution and high throughput longitudinal recording systems allowed the acquisition of faster and more detailed images at mesoscale level for months on large cohorts of animals. This generates large amounts of data that need to be stored, pre-processed and analyzed in a timely manner which can become a limiting factor as experiments become larger and more complex.   

The ***umIToolbox***, or ***umIT*** for short is a toolbox written in MATLAB created to help researchers manage and automate critical steps of processing large-scale imaging datasets. The toolbox is currently in active development and new features are added periodically. The toolbox can be found in our [GitHub repository](https://github.com/S-Belanger/Umit). If you encounter any bugs or have suggestions on how to improve the toolbox, let us know by opening an issue in the [GitHub issue tracker](https://github.com/S-Belanger/Umit/issues). You can also use the project's [discussion forum](https://github.com/S-Belanger/Umit/discussions) to share ideas and questions with the developers and other users of the toolbox.

This toolbox is dedicated to our beloved friend ***Umit Keysan*** (1992-2019).

___

### Getting started with ***umIT***
Here, you will find the minimal information to get started with the toolbox.
* [An overview](#an-overview)
* [umIT installation](/umit_install.md)
* [The main user interface](/mainGUI.md)
* [Creating a new project](/ht_create_new_project.md) (in construction)
* [Running an analysis pipeline](/ht_run_pipeline.md) (in construction)
* [Visualizing results](/ht_viz_data.md)  (in construction)

___

### User Guides
* [Build a data processing pipeline](/ht_buildpipeline.md)
* [Extract data from ROI](/ht_getdatafromroi.md) using the standalone version of [*DataViewer*](/dataviewer.md) (link not functional)
* [Import Data from LabeoTech Systems](/docs/userDocs/dataimportlabeo.md)

___

### Documentation

Here you will find more information on the features of each **app** from the toolbox and the documentation of the built-in **analysis functions**.   

Some of the apps work as standalone, meaning that they can operate without using the main interface (i.e., the experiment management system). This is particularly useful if you want to quickly analyse and visualize a small number of imaging recording sessions before running it through the experiment management system.
##### Apps
* [DataViewer](/dataviewer.md) : Basic exploration of imaging data. Works as standalone.
* [ROImanager](/roimanager.md) : Creation and management of Regions of interest. Works as standalone.

##### Analysis functions
* [run_ImagesClassification](/docs/devDocs/run_imagesclassification.md)

##### *umIT*'s built-in functions and classes
* [loadDatFile](/docs/devDocs/loaddatfile.md)

___
## An overview
A typical imaging **project** consists of one or more cohorts of **subject**s (e.g. mice) that undergo one or more **acquisition** (i.e. recording) sessions. Frequently, other recording **modalities** are associated with the imaging data such as behavioral responses, eye/body tracking, etc. The toolbox follows the same organization principle where one can manage subjects, acquisitions and recording modalities for a given project.   

![fig1](/assets/img/umIT_concept_org_img_exp.png)
###### Fig 1. Conceptual organization of a typical imaging project.

The management of the project datasets is done through a set of *objects* that control the project (class Protocol), subjects (class Subject), acquisitions (class Acquisition) and modalities (abstract class modality and sub classes such as *FluorescenceImaging* class, for imaging data). Once the project file is created, the toolbox can automatically detect new subjects or acquisitions and update the file.   

In addition to the experiment management tool, the toolbox provides a series of **analysis functions** and a **pipeline** engine that allows the automation of several steps of the data processing. These two modules (experiment management and data analysis) were created to be adaptable and easily extensible in order to fulfill the needs of different neuroimaging experiment designs. Finally, ***umIT*** provides tools to visualize imaging data and to perform statistical comparisons between groups.

![fig2](/assets/img/umIT_4axis.png)
###### Fig 2. The 4 main roles of **umIT**.  
